/**
 * @file       inflater.cpp
 * @date       Jun 2, 2020
 * @author     Martin Rizzo | <martinrizzo@gmail.com>
 * @copyright  Copyright (c) 2020 Martin Rizzo.
 *             This project is released under the MIT License.
 * -------------------------------------------------------------------------
 *  Inflater - One-header library to decode data compressed with the Deflate algorithm.
 * -------------------------------------------------------------------------
 *  Copyright (c) 2020 Martin Rizzo
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 *  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 *  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * -------------------------------------------------------------------------
 */
#include "inflater.h"
#include <stdlib.h>
#include <string.h>
#include "PDZip_ReversedHuffmanDecoder.h"

#include <stdio.h>

static size_t min(size_t a, size_t b) { return a<b ? a : b; }

#define EndOfBlock           256
#define MaxValidLengthCode   285
#define MaxValidDistanceCode 29


//======================================================================================================================
#   pragma mark - LOW-LEVEL BIT STREAM

struct BitStream {
  
    BitStream(const unsigned char* inputPtr, const unsigned char* inputEnd, unsigned bitbuffer, unsigned bitcounter) :
    _inputPtr(inputPtr),
    _inputEnd(inputEnd),
    _bitbuffer(bitbuffer),
    _bitcounter(bitcounter)
    {
        assert( inputPtr!=NULL );
        assert( inputEnd!=NULL );
        assert( 0<=bitcounter && bitcounter<=32 );
    }
    
    inline bool loadNextByte() {
        if (_inputPtr==_inputEnd) { return false; }
        _bitbuffer  |= (*_inputPtr++ << _bitcounter);
        _bitcounter += 8;
        return true;
    }
    
    inline bool readBits(unsigned* dest, int numberOfBitsToRead) {
        assert( dest!=NULL );
        while (_bitcounter<numberOfBitsToRead) { if (!loadNextByte()) { return false; }  }
        (*dest) = _bitbuffer & ((1<<numberOfBitsToRead)-1);
        _bitbuffer  >>= numberOfBitsToRead;
        _bitcounter  -= numberOfBitsToRead;
        return true;
    }
    
    inline bool readBits(unsigned* dest, const PDZip::ReversedHuffmanDecoder* decoder) {
        assert( dest!=NULL );
        assert( decoder!=NULL );
        assert( decoder->isLoaded() );
        PDZip::ReversedHuffmanDecoder::Data data;

        // TODO: optimize `readBits` for when when there are enough input buffer data
        if ( _bitcounter>0 ) {
            data = decoder->decode8(_bitbuffer);
            if (data.isValid() && data.length()<=_bitcounter) {
                _bitbuffer >>= data.length();
                _bitcounter -= data.length();
                (*dest) = data.code();
                return true;
            }
        }
        if ( _bitcounter<8 && !loadNextByte() ) { return false; }
        data = decoder->decode8(_bitbuffer);
        if (data.isValid()) {
            _bitbuffer  >>= data.length();
            _bitcounter  -= data.length();
            (*dest) = data.code();
            return true;
        }
        if ( _bitcounter<16 && !loadNextByte() ) { return false; }
        data = decoder->decode16(_bitbuffer,data);
        if (data.isValid()) {
            _bitbuffer >>= data.length();
            _bitcounter -= data.length();
            (*dest) = data.code();
            return true;
        }
        assert( false );
        (*dest) = static_cast<unsigned>(-1);
        return true;
    }
    
    inline bool readWord(unsigned* dest) {
        static const unsigned WordBits = 16; // < number of bits in a Word
        assert( dest!=NULL );
        assert( 8*sizeof(*dest)>=WordBits );
        // skip any padding bits (bitbuffer must be byte aligned before reading a Word)
        const unsigned  numberOfPaddingBits = (_bitcounter%8);
        assert( numberOfPaddingBits<=_bitcounter );
        _bitbuffer  >>= numberOfPaddingBits;
        _bitcounter  -= numberOfPaddingBits;
        while ( _bitcounter<WordBits ) { if (!loadNextByte()) { return false; } }
        assert( _bitcounter==WordBits );
        (*dest) = (_bitbuffer>>8 & 0x00FF) | (_bitbuffer<<8 & 0xFF00);
        _bitbuffer  = 0;
        _bitcounter = 0;
        return true;
    }
    
    inline bool readDWord(unsigned* dest) {
        static const unsigned DWordBits = 32; // < number of bits in a DWord
        assert( dest!=NULL );
        assert( 8*sizeof(*dest)>=DWordBits );
        // skip any padding bits (bitbuffer must be byte aligned before reading a Word)
        const unsigned  numberOfPaddingBits = (_bitcounter%8);
        _bitbuffer  >>= numberOfPaddingBits;
        _bitcounter  -= numberOfPaddingBits;
        while ( _bitcounter<DWordBits ) { if (!loadNextByte()) { return false; } }
        assert( _bitcounter==DWordBits );
        (*dest) = (_bitbuffer>>24 & 0x000000FF) | (_bitbuffer>>8 & 0x0000FF00) | (_bitbuffer<<8 & 0x00FF0000) | (_bitbuffer<<24 & 0xFF000000);
        _bitbuffer  = 0;
        _bitcounter = 0;
        return true;
    }

    inline bool readBytes(unsigned char* outputPtr, size_t* inout_numberOfBytes) {
        assert( outputPtr!=NULL );
        assert( inout_numberOfBytes!=NULL );
        assert( _bitcounter==0 );
        
        const size_t numberOfBytesToRead = (*inout_numberOfBytes);
        const size_t successfulBytes     = min( (_inputEnd-_inputPtr), numberOfBytesToRead );
        memcpy(outputPtr, _inputPtr, successfulBytes);
        _inputPtr += successfulBytes;
        (*inout_numberOfBytes) = successfulBytes;
        return (successfulBytes == numberOfBytesToRead);
    }
    
    unsigned             _bitbuffer;
    unsigned             _bitcounter;
    const unsigned char* _inputPtr;
    const unsigned char* _inputEnd;
};


//======================================================================================================================
#   pragma mark - ???


typedef enum InfMode {
    InfMode_Uninitialized,
    InfMode_Take,
    InfMode_Feed
} InflaterMode;

typedef enum InfFlags {
    InfFlags_FreeOutputBuffer,
    InfFlags_FreeInputBuffer,
    InfFlags_FreeItself
} InfFlags;

typedef enum BlockType {
    BlockType_Uncompressed = 0,
    BlockType_FixedHuffman = 1,
    BlockType_DynamicHuffman = 2,
} BlockType;


 struct CodeLengthsReader {
     unsigned           _command;
     unsigned           _code;
     unsigned           _length;
     unsigned           _repetitions;
     PDZip::CodeLengths _codeLengths;

     void reset() { _command = _length = _repetitions = 0; }
     void begin() { _code=0; _codeLengths.open(); }
     void end()   { _codeLengths.close(); }
     const PDZip::CodeLengths& operator *() const { return _codeLengths; }

     bool read(void* bstream_ptr, unsigned numberOfCodes) {
         assert( _codeLengths.isOpen() );
         assert( numberOfCodes>0 );
         assert( bstream_ptr!=NULL );
         static const unsigned MaximumNumberOfCodes = 19;
         static const unsigned order[MaximumNumberOfCodes] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };
         
         BitStream&           bstream       = *static_cast<BitStream*>(bstream_ptr);
         unsigned char* const internalArray = _codeLengths.getInternalArray(MaximumNumberOfCodes);
         
         const unsigned end = (numberOfCodes<MaximumNumberOfCodes) ? numberOfCodes : MaximumNumberOfCodes;
         while (_code<end) {
             unsigned length;
             if ( !bstream.readBits(&length,3) ) { return false; }
             internalArray[order[_code]] = length;
             ++_code;
         }
         while (_code<MaximumNumberOfCodes) { internalArray[order[_code++]] = 0; }
         _codeLengths.addFromInternalArray(internalArray,MaximumNumberOfCodes);
         return true;
     }
     
     bool read(void* bstream_ptr, unsigned numberOfCodes, const PDZip::ReversedHuffmanDecoder* decoder) {
         assert( _codeLengths.isOpen() );
         assert( numberOfCodes>0 );
         assert( bstream_ptr!=NULL );
         static const unsigned MaxValidLength = 15; // < the maximum length that can be read
         enum Command {
             Command_ReadNext            = 0,  // < load next command
             Command_CopyPreviousLength  = 16, // < repeat the previous length
             Command_RepeatZeroLength_3  = 17, // < repeat zero length (3 or more times)
             Command_RepeatZeroLength_11 = 18  // < repeat zero length (11 or more times)
         };
         BitStream& bstream = *static_cast<BitStream*>(bstream_ptr);
         unsigned value = 0;

         // IMPORTANT: add any repetition that is pending from a previous load (ex: literals-table > distance-table)
         while (_code<numberOfCodes && _repetitions>0) {
             _codeLengths.add(_code++,_length); --_repetitions;
         }
         
         // add (one by one) all codes with their codeLengths taking into account the number of repetitions
         while ( _code<numberOfCodes ) {
             assert( _repetitions==0 );
             
             // read a new command
             if ( _command==Command_ReadNext ) {
                 if ( !bstream.readBits(&_command, decoder) ) { return false; }
             }
             // process simple command
             if ( _command<=MaxValidLength ) {
                 _codeLengths.add(_code++,_length=_command);
             }
             // process command with repetition
             else {
                 switch (_command) {
                     case Command_CopyPreviousLength:
                         if (_code==0) { /* codeLengths.markAsError(); */ return true; }
                         if ( !bstream.readBits(&value,2) ) { return false; }
                         _repetitions = 3 + value;
                         break;
                     case Command_RepeatZeroLength_3:
                         if ( !bstream.readBits(&value,3) ) { return false; }
                         _length      = 0;
                         _repetitions = 3 + value;
                         break;
                     case Command_RepeatZeroLength_11:
                         if ( !bstream.readBits(&value,7) ) { return false; }
                         _length      = 0;
                         _repetitions = 11 + value;
                         break;
                     default:
                         /* codeLengths.markAsError(); */
                         return true;
                 }
                 while (_code<numberOfCodes && _repetitions>0) {
                     _codeLengths.add(_code++,_length); --_repetitions;
                 }
             }
             // mark current command as finished
             _command = Command_ReadNext;
         }
         return true;

     }
 };


#   define op_goto(dest_step)        step=dest_step; break;
#   define op_fallthrough(next_step) step=next_step;
#   define op_return(x)                                       exit_status=x; break;
#   define op_end(x)                 step=InfStep_End;        exit_status=x; break;
#   define op_fatal_error(x)         step=InfStep_FatalError; exit_status=(InfAction)x; break;
    
typedef enum InfStep {
    InfStep_Start,
    InfStep_ReadZLibHeader,
    InfStep_ReadZLibHeader2,
    InfStep_ProcessNextBlock,
    InfStep_ReadBlockHeader,
    
    InfStep_Load_FixedHuffmanDecoders,
    InfStep_Load_DynamicHuffmanDecoders,
    InfStep_Load_FrontHuffmanTable,
    InfStep_Load_FrontHuffmanTable2,
    InfStep_Load_LiteralHuffmanTable,
    InfStep_Load_LiteralHuffmanTable2,
    InfStep_Load_DistanceHuffmanTable,
    InfStep_Load_DistanceHuffmanTable2,
    
    InfStep_Process_UncompressedBlock,
    InfStep_Process_CompressedBlock,
    
    InfStep_Read_LiteralOrLength,
    InfStep_Read_LiteralOrLength2,
    InfStep_Read_LengthBits,
    InfStep_Read_Distance,
    InfStep_Read_DistanceBits,
    
    InfStep_Output_UncompressedBlock,
    InfStep_Output_RepeatedSequence,
    
    InfStep_FatalError,
    InfStep_End
} InfStep;


static PDZip::ReversedHuffmanDecoder* getFixedLiteralDecoder() {
    static PDZip::ReversedHuffmanDecoder decoder;
    if ( !decoder.isLoaded() ) {
        PDZip::CodeLengths codeLengths;
        codeLengths.open();
        codeLengths.addRange(  0,144, 8);
        codeLengths.addRange(144,256, 9);
        codeLengths.addRange(256,280, 7);
        codeLengths.addRange(280,288, 8);
        codeLengths.close();
        decoder.load(codeLengths);
    }
    return &decoder;
}

static PDZip::ReversedHuffmanDecoder* getFixedDistanceDecoder() {
    static PDZip::ReversedHuffmanDecoder decoder;
    if ( !decoder.isLoaded() ) {
        PDZip::CodeLengths codeLengths;
        codeLengths.open();
        codeLengths.addRange(0,32, 5);
        codeLengths.close();
        decoder.load(codeLengths);
    }
    return &decoder;
}

/* TODO: remove these globals! */
static CodeLengthsReader              _codeLengths;
static PDZip::ReversedHuffmanDecoder  _huffmanDecoder0;
static PDZip::ReversedHuffmanDecoder  _huffmanDecoder1;
static PDZip::ReversedHuffmanDecoder  _huffmanDecoder2;


/**
 * Decompress a chunk of compressed data
 *
 * @param inputPtr
 *     The pointer to the position from where the compressed data will be read
 *
 * @param inout_inputBytes
 *     [in]  The number of bytes available to read.
 *     [out] The number of bytes successfully read
 *
 * @param outputPtr
 *     The pointer to the position where the resulting uncompressed data will be written to.
 *     This pointer must point to a valid position within the output buffer.
 *
 * @param inout_outputBytes
 *     [in]  The number of bytes avariable to write.
 *     [out] The number of bytes successfully written
 *
 * @param outputBufferBegin
 *     The beginning of the output buffer (used when copying previous sequences)
 *
 * @param outputBufferSize
 *     The total size of the output buffer (in number of bytes)
 */
InfAction Inf_decompress(Inflater*                  inf,
                         const unsigned char* const inputPtr,
                         size_t*                    inout_inputBytes,
                         unsigned char* const       outputPtr,
                         size_t*                    inout_outputBytes,
                         unsigned char* const       outputBufferBegin,
                         size_t                     outputBufferSize)
{
    assert( inputPtr!=NULL );
    assert( inout_inputBytes!=NULL && (*inout_inputBytes)>0 );
    assert( outputPtr!=NULL );
    assert( inout_outputBytes!=NULL && (*inout_outputBytes)>0 );
    assert( outputBufferBegin!=NULL );
    assert( outputBufferSize>=(32*1024) );

    static const unsigned int lengthStarts[] = {
         3,  4,  5,  6,  7,  8,  9,  10,  11,  13,  15,  17,  19,  23,  27,  31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
    static const unsigned int lengthExtraBits[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
    static const unsigned int distanceStarts[] = {
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577};
    static const unsigned int distanceExtraBits[] = {
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
        12, 12, 13, 13};

    BitStream bstream(inputPtr, inputPtr+(*inout_inputBytes), inf->_bitbuffer, inf->_bitcounter);
    unsigned char* const writeEnd = outputPtr + (*inout_outputBytes);
    unsigned char*       writePtr = outputPtr;

    unsigned       temp;
    size_t         numberOfBytes;
    bool           canRead;
    unsigned char* sequencePtr;
    
    
    InfStep   step        = (InfStep)inf->_op;
    InfAction exit_status = (InfAction)(step!=InfStep_End && step!=InfStep_FatalError ? InfAction_ProcessNextChunk : InfError_Failed);
    while ( exit_status==InfAction_ProcessNextChunk ) {
        switch (step) {
                
            case InfStep_Start:
                inf->_lastBlock = false;
                op_fallthrough(InfStep_ReadZLibHeader);
                
            //-------------------------------------------------------------------------------------
            // Read_ZLibHeader
            //
            case InfStep_ReadZLibHeader:
                if ( !bstream.readWord(&inf->_zlibheader) ) { op_return(InfAction_FillInputBuffer); }
                op_fallthrough(InfStep_ReadZLibHeader2);
                
            case InfStep_ReadZLibHeader2:
                if ( (inf->_zlibheader % 31)!=0 ) { op_fatal_error(InfError_BadZLibHeader); }
                inf->_zlibheader_method = inf->_zlibheader>>8 & 0x0F;
                inf->_zlibheader_wsize  = 1 << (8 + (inf->_zlibheader>>12 & 0x0F));
                inf->_zlibheader_level  = inf->_zlibheader>>6 & 0x03;
                if (inf->_zlibheader_method!=8)       { op_fatal_error(InfError_UnsupportedZLibHeader); }
                if (inf->_zlibheader&0x20)            { op_fatal_error(InfError_UnsupportedZLibHeader); }
                if (inf->_zlibheader_wsize>(32*1024)) { op_fatal_error(InfError_UnsupportedZLibHeader); }
                //if (_zlibheader_wsize>wrappingBufferSize) { op_fatal_error(Status_UnsupportedZLibHeader); }
                op_fallthrough(InfStep_ProcessNextBlock);
                
            //-------------------------------------------------------------------------------------
            // Read_BlockHeader
            //
            case InfStep_ProcessNextBlock:
                if (inf->_lastBlock) {
                    printf(" > END OF STREAM\n\n");
                    op_end(InfAction_Finish);
                }
                op_fallthrough(InfStep_ReadBlockHeader);
                
            case InfStep_ReadBlockHeader:
                if ( !bstream.readBits(&temp,3) ) { op_return(InfAction_FillInputBuffer); }
                inf->_lastBlock = (temp&0x01)==0x01;
                inf->_blocktype = (temp>>1);
                if      (inf->_blocktype==BlockType_Uncompressed)   {
                    printf(" > Start Uncompressed Block\n");
                    op_goto(InfStep_Process_UncompressedBlock);
                }
                else if (inf->_blocktype==BlockType_FixedHuffman)   {
                    printf(" > Start Fixed Huffman Block\n");
                    op_goto(InfStep_Load_FixedHuffmanDecoders);
                }
                else if (inf->_blocktype==BlockType_DynamicHuffman) {
                    printf(" > Start Dynamic Huffman Block\n");
                    op_goto(InfStep_Load_DynamicHuffmanDecoders);
                }
                printf(" > Fatal Error (Bad Block Type)\n");
                op_fatal_error(InfError_BadBlockType);
                
            //-------------------------------------------------------------------------------------
            // Process_UncompressedBlock
            //
            case InfStep_Process_UncompressedBlock:
                if ( !bstream.readDWord(&temp) ) { op_return(InfAction_FillInputBuffer); }
                inf->_seq_len = (temp >> 16);
                if ( inf->_seq_len != ((~temp)&0xFFFF) ) { op_fatal_error(InfError_BadBlockLength); }
                op_fallthrough(InfStep_Output_UncompressedBlock);
                
            case InfStep_Output_UncompressedBlock:
                numberOfBytes = min( (writeEnd-writePtr), inf->_seq_len );
                canRead       = bstream.readBytes(writePtr,&numberOfBytes);
                inf->_seq_len -= numberOfBytes;
                writePtr      += numberOfBytes;
                if ( inf->_seq_len>0 ) { op_return(InfAction_UseOutputBufferContent); }
                else if ( !canRead )   { op_return(InfAction_FillInputBuffer);        }
                
                assert( inf->_seq_len==0 );
                op_goto(InfStep_ProcessNextBlock);

            //-------------------------------------------------------------------------------------
            // Load FixedHuffmanDecoders
            //
            case InfStep_Load_FixedHuffmanDecoders:
                inf->_literalDecoder  = getFixedLiteralDecoder();
                inf->_distanceDecoder = getFixedDistanceDecoder();
                op_goto(InfStep_Process_CompressedBlock);
                
            //-------------------------------------------------------------------------------------
            // Load DynamicHuffmanDecoders
            //
            case InfStep_Load_DynamicHuffmanDecoders:
                if ( !bstream.readBits(&temp,5+5+4) ) { op_return(InfAction_FillInputBuffer); }
                inf->_hlit  = temp     & 0x1F;
                inf->_hdist = temp>>5  & 0x1F;
                inf->_hclen = temp>>10 & 0x0F;
                op_fallthrough(InfStep_Load_FrontHuffmanTable);

                //----------------------------------
                // load "front" huffman table
                //
            case InfStep_Load_FrontHuffmanTable:
                _codeLengths.reset();
                _codeLengths.begin();
                op_fallthrough(InfStep_Load_FrontHuffmanTable2);
                
            case InfStep_Load_FrontHuffmanTable2:
                if ( !_codeLengths.read(&bstream,inf->_hclen+4) ) { op_return(InfAction_FillInputBuffer); }
                _codeLengths.end();
                inf->_frontDecoder = &_huffmanDecoder0;
                inf->_frontDecoder->load( *_codeLengths );
                op_fallthrough(InfStep_Load_LiteralHuffmanTable);
                
                //----------------------------------
                // load literal-length huffman table
                //
            case InfStep_Load_LiteralHuffmanTable:
                _codeLengths.reset();
                _codeLengths.begin();
                op_fallthrough(InfStep_Load_LiteralHuffmanTable2);
                
            case InfStep_Load_LiteralHuffmanTable2:
                if ( !_codeLengths.read(&bstream,inf->_hlit+257,inf->_frontDecoder) ) { op_return(InfAction_FillInputBuffer); }
                _codeLengths.end();
                inf->_literalDecoder = &_huffmanDecoder1;
                inf->_literalDecoder->load( *_codeLengths );
                op_fallthrough(InfStep_Load_DistanceHuffmanTable);
                
                //----------------------------------
                // load distance huffman table
                //
            case InfStep_Load_DistanceHuffmanTable:
                _codeLengths.begin();
                op_fallthrough(InfStep_Load_DistanceHuffmanTable2);
                
            case InfStep_Load_DistanceHuffmanTable2:
                if ( !_codeLengths.read(&bstream,inf->_hdist+1,inf->_frontDecoder) ) { op_return(InfAction_FillInputBuffer); }
                _codeLengths.end();
                inf->_distanceDecoder = &_huffmanDecoder2;
                inf->_distanceDecoder->load( *_codeLengths );
                op_fallthrough(InfStep_Process_CompressedBlock);
                
            //-------------------------------------------------------------------------------------
            // Process Compressed Block
            //
            case InfStep_Process_CompressedBlock:
            case InfStep_Read_LiteralOrLength:
                if ( !bstream.readBits(&inf->_literal,inf->_literalDecoder) ) { op_return(InfAction_FillInputBuffer); }
                op_fallthrough(InfStep_Read_LiteralOrLength2);
                
            case InfStep_Read_LiteralOrLength2:
                if (inf->_literal <EndOfBlock) {
                    if (writePtr==writeEnd) { op_return(InfAction_UseOutputBufferContent); }
                    *writePtr++ = inf->_literal;
                    op_goto(InfStep_Read_LiteralOrLength);
                }
                else if (inf->_literal==EndOfBlock) {
                    printf(" > EndOfBlock\n");
                    op_goto(InfStep_ProcessNextBlock);
                }
                else if (inf->_literal>MaxValidLengthCode) { op_fatal_error(InfError_BadBlockContent); }
                inf->_literal -= 257;
                op_fallthrough(InfStep_Read_LengthBits);
                
            case InfStep_Read_LengthBits:
                if ( !bstream.readBits(&temp, lengthExtraBits[inf->_literal]) ) { op_return(InfAction_FillInputBuffer); }
                inf->_seq_len   =   temp + lengthStarts[inf->_literal];
                assert( inf->_seq_len<(32*1024) );
                op_fallthrough(InfStep_Read_Distance);
                
            case InfStep_Read_Distance:
                if ( !bstream.readBits(&inf->_literal,inf->_distanceDecoder) ) { op_return(InfAction_FillInputBuffer); }
                if (inf->_literal>MaxValidDistanceCode) { op_fatal_error(InfError_BadBlockContent); }
                op_fallthrough(InfStep_Read_DistanceBits);
                
            case InfStep_Read_DistanceBits:
                if ( !bstream.readBits(&temp, distanceExtraBits[inf->_literal]) ) { op_return(InfAction_FillInputBuffer); }
                inf->_seq_dist =  temp + distanceStarts[inf->_literal];
                assert( inf->_seq_dist<(32*1024) );
                op_fallthrough(InfStep_Output_RepeatedSequence);

            //-------------------------------------------------------------------------------------
            // Output_RepeatedSequence
            //
            //         outputBegin     writePtr                      writeEnd
            //              V             V                             V
            // <... ghost > # [[ straight * overlapped  ...... ghost ]] #
            //
            case InfStep_Output_RepeatedSequence:
                sequencePtr = writePtr - inf->_seq_dist;
                if (sequencePtr<outputBufferBegin) {
                    sequencePtr += outputBufferSize;
                    //-- copy bytes ---
                    inf->_seq_len -= numberOfBytes = min( inf->_seq_len, (writeEnd-sequencePtr) );
                    memmove( writePtr, sequencePtr, numberOfBytes ); writePtr+=numberOfBytes;
                    //-----------------
                    if ( inf->_seq_len==0 ) { op_goto(InfStep_Read_LiteralOrLength); }
                    if ( writePtr==writeEnd )  { op_return(InfAction_UseOutputBufferContent); }
                    sequencePtr = outputBufferBegin;
                }
                //-- copy bytes ---
                inf->_seq_len -= numberOfBytes = min( inf->_seq_len, (writeEnd-writePtr) );
                while ( numberOfBytes-->0 ) { *writePtr++ = *sequencePtr++; }
                //-----------------
                if ( inf->_seq_len==0 ) { op_goto(InfStep_Read_LiteralOrLength); }
                if ( writePtr==writeEnd )  { op_return(InfAction_UseOutputBufferContent); }
                op_fatal_error(InfError_BadBlockContent);
                
                
            case InfStep_FatalError:
            case InfStep_End:
                break;
        }
    }
    
    inf->_op         = step;
    inf->_bitbuffer  = bstream._bitbuffer;
    inf->_bitcounter = bstream._bitcounter;
    
    printf(" # exit_stats = %d\n", exit_status);

    (*inout_inputBytes)  = (bstream._inputPtr - inputPtr);
    (*inout_outputBytes) = (writePtr          - outputPtr);
    return exit_status;
}





Inflater* inflaterCreate(void* workingBuffer, size_t workingBufferSize) {
    Inflater* inf = (Inflater*)malloc(sizeof(Inflater));
    
    /*
    inf->obj = new PDZip::Inflater;
    inf->obj->init();
    */
    
    inf->mode            = InfMode_Uninitialized;
    inf->flags           = InfFlags_FreeItself;
    inf->action          = InfAction_Init;
    inf->error           = InfError_None;
    
    inf->finished = 0;
    
    inf->userPtr          = NULL;
    inf->dataProviderFunc = NULL;
    inf->dataReceiverFunc = NULL;
    
    inf->inputChunk      = NULL;
    inf->inputChunkSize  = 0;
    inf->outputChunk     = NULL;
    inf->outputChunkSize = 0;
    inf->outputBufferContentSize = 0;
    
    /*---- inflaterTake / inflaterFeed ----------- */
    inf->helperOutput.buffer     = NULL;
    inf->helperOutput.bufferSize = 0;
    inf->helperInput.buffer      = NULL;
    inf->helperInput.bufferSize  = 0;
    inf->takeOutputPtr           = NULL;
    inf->takeOutputRemaining     = 0;
    inf->providedData.buffer     = NULL;
    inf->providedData.bufferSize = 0;
    
    /*---- decompress ---------------------------- */
    inf->_op         = 0;
    inf->_bitbuffer  = 0;
    inf->_bitcounter = 0;
    
    
    if (workingBuffer!=NULL && workingBufferSize>=InfHelperBufferSize) {
        inf->helperOutput.buffer     = workingBuffer;
        inf->helperOutput.bufferSize = workingBufferSize;
    }
    return inf;
}


/**
 * Decompresses a chunk of data
 * @param inflater          Pointer to the `Inflater` object created with `inflaterCreate(..)`
 * @param outputBuffer      Pointer to the destination buffer where decompressed data will be stored
 * @param outputBufferSize  The available capacity of `destBuffer` in number of bytes
 * @param inputBuffer       Pointer to the source buffer from where compressed data is read
 * @param inputBufferSize   The length of `sourBuffer` in number of bytes
 * @returns
 *     current status
 */
InfAction inflaterProcessChunk(Inflater*   inflater,
                               void*       outputBuffer,
                               size_t      outputBufferSize,
                               const void* inputBuffer,
                               size_t      inputBufferSize) {

    InfAction status;
    size_t inputBytes, outputBytes;
    const Byte* const inputBufferEnd  = (Byte*)inputBuffer  + inputBufferSize;
    Byte* const       outputBufferEnd = (Byte*)outputBuffer + outputBufferSize;
    
    assert( inputBuffer !=NULL && inputBufferSize >0 );
    assert( outputBuffer!=NULL && outputBufferSize>0 );
    
    if (inflater->finished) {
        inflater->action=InfAction_Finish;
    }
    
    switch ( inflater->action ) {
            
        case InfAction_Init:
            inflater->inputChunk  = (const Byte*)inputBuffer;
            inflater->outputChunk = (Byte      *)outputBuffer;
            break;
            
        case InfAction_FillInputBuffer:
            inflater->inputChunk   = (const Byte*)inputBuffer;
            inflater->outputChunk += inflater->outputChunkSize;
            break;
            
        case InfAction_UseOutputBufferContent:
            inflater->inputChunk  += inflater->inputChunkSize;
            inflater->outputChunk  = (Byte*)outputBuffer;
            inflater->outputBufferContentSize = 0;
            break;
            
        case InfAction_ProcessNextChunk:
            inflater->inputChunk  += inflater->inputChunkSize;
            inflater->outputChunk += inflater->outputChunkSize;
            break;
            
        /* InfAction_Finish */
        default:
            inflater->inputChunkSize  = 0;
            inflater->outputChunkSize = 0;
            return inflater->action;
    }
    
    assert( inflater->inputChunk  < inputBufferEnd  );
    assert( inflater->outputChunk < outputBufferEnd );
    
    inputBytes  = (inputBufferEnd  - inflater->inputChunk );
    outputBytes = (outputBufferEnd - inflater->outputChunk);
    
    status = Inf_decompress(inflater,
                            inflater->inputChunk , &inputBytes ,
                            inflater->outputChunk, &outputBytes,
                            (unsigned char*)outputBuffer, outputBufferSize);
    
    /*
    status = inflater->obj->decompress(inflater->inputChunk , &inputBytes ,
                                       inflater->outputChunk, &outputBytes,
                                       (unsigned char*)outputBuffer, outputBufferSize);
    */
    
    inflater->action          = status>=0 ? (InfAction)status : InfAction_Finish;
    inflater->error           = status <0 ? (InfError)status  : InfError_None;
    inflater->inputChunkSize  = inputBytes;
    inflater->outputChunkSize = outputBytes;
    inflater->outputBufferContentSize += outputBytes;
    
    /* hack to add an 'UseOutputBufferContent' before 'Finish' */
    inflater->finished = (inflater->action==InfAction_Finish);
    if (inflater->finished && outputBytes>0) {
        inflater->action=InfAction_UseOutputBufferContent;
    }
    
    return inflater->action;
}


#define isPaused(inf) (inf->action==InfAction_FillInputBuffer && inf->providedData.bufferSize==0)
#define isFinished(inf) (inf->action==InfAction_Finish)

/**
 * Decompresses a chunk of data
 * @param inf                   Pointer to the `Inflater` object created with `inflaterCreate(..)`
 * @param decompressedData      Pointer to the destination buffer where decompressed data will be stored
 * @param decompressedDataSize  The available capacity of `destBuffer` in number of bytes
 */
size_t inflaterTake(Inflater* inf, void* decompressedData, size_t decompressedDataSize) {
    /* InfAction action; */
    size_t sourAvailableBytes,destBytesStillNeeded;
    Byte* dest; const Byte* sour;
    Byte* outputBuffer; size_t outputBufferSize;
    int fillInputBuffer;

    /* initialization */
    assert( (inf->mode==InfMode_Uninitialized) || (inf->mode==InfMode_Take) );
    if (inf->mode==InfMode_Uninitialized) {
        inf->mode= InfMode_Take;
        if (inf->helperOutput.buffer==NULL) {
            inf->helperOutput.buffer = malloc( inf->helperOutput.bufferSize=InfHelperBufferSize );
            inf->flags |= InfFlags_FreeOutputBuffer;
        }
        if (inf->helperInput.buffer==NULL) {
            inf->helperInput.buffer  = malloc( inf->helperInput.bufferSize=InfHelperBufferSize );
            inf->flags |= InfFlags_FreeInputBuffer;
        }
        inf->takeOutputPtr           = (Byte*)inf->helperOutput.buffer;
        inf->takeOutputRemaining     = 0;
        inf->providedData.buffer     = NULL;
        inf->providedData.bufferSize = 0;
    }

    outputBuffer         = (Byte*)inf->helperOutput.buffer;
    outputBufferSize     =        inf->helperOutput.bufferSize;
    sour                 =        inf->takeOutputPtr;
    sourAvailableBytes   =        inf->takeOutputRemaining;
    dest                 = (Byte*)decompressedData;
    destBytesStillNeeded =        decompressedDataSize;
    
    
    if (isFinished(inf)) { return 0; }

    while ( destBytesStillNeeded>0 )
    {
        fillInputBuffer = (inf->action==InfAction_Init || inf->action==InfAction_FillInputBuffer);
        
        if (inf->action==InfAction_UseOutputBufferContent && sourAvailableBytes>0 ) {
            /* take available bytes */
            const size_t size = sourAvailableBytes<destBytesStillNeeded ? sourAvailableBytes : destBytesStillNeeded;
            memcpy(dest,sour,size); dest+=size; destBytesStillNeeded-=size; sour+=size; sourAvailableBytes-=size;
            if (sourAvailableBytes==0) { fillInputBuffer=1; sour=outputBuffer; }
        }
        if ( fillInputBuffer ) {
            /* get compressed data from the data provider (when required) */
            if (inf->action==InfAction_Init || inf->action==InfAction_FillInputBuffer) {
                inf->dataProviderFunc(inf, &(inf->providedData=inf->helperInput) );
                if ( inf->providedData.bufferSize==0 ) { break; }
            }
            /* decompress the provided data and generate more sourAvailableBytes */
            do {
                inflaterProcessChunk(inf, outputBuffer, outputBufferSize, inf->providedData.buffer, inf->providedData.bufferSize);
                sourAvailableBytes += inf->outputChunkSize;
            } while ( inf->action==InfAction_ProcessNextChunk );
        }
    }
    inf->takeOutputPtr=sour; inf->takeOutputRemaining=sourAvailableBytes;
    return (decompressedDataSize-destBytesStillNeeded);
}

/**
 * Decompresses a chunk of data
 * @param inf                 Pointer to the `Inflater` object created with `infalterCreate(..)`
 * @param compressedData      Pointer to the source buffer from where compressed data is read
 * @param compressedDataSize  The length of `compressedData` in number of bytes
 */
size_t inflaterFeed(Inflater* inf, const void* compressedData, size_t compressedDataSize) {
    InfAction action; size_t numberOfConsumedBytes=0;
    Byte* outputBuffer; size_t outputBufferSize;
    assert( (inf->mode==InfMode_Uninitialized) || (inf->mode=InfMode_Feed) );
    assert( inf->dataReceiverFunc!=NULL );

    /* initialization */
    if (inf->mode==InfMode_Uninitialized) {
        inf->mode= InfMode_Feed;
        if (inf->helperOutput.buffer==NULL) {
            inf->helperOutput.buffer = (Byte*)malloc( inf->helperOutput.bufferSize=InfHelperBufferSize );
            inf->flags = InfFlags_FreeOutputBuffer;
        }
    }
    
    /* decompression */
    outputBuffer     = (Byte*)inf->helperOutput.buffer;
    outputBufferSize =        inf->helperOutput.bufferSize;
    do {
        action = inflaterProcessChunk(inf, outputBuffer, outputBufferSize, compressedData, compressedDataSize);
        if (action==InfAction_UseOutputBufferContent) {
            inf->dataReceiverFunc(inf, outputBuffer, inf->outputBufferContentSize);
        }
        numberOfConsumedBytes += inf->inputChunkSize;
    } while (action==InfAction_ProcessNextChunk || action==InfAction_UseOutputBufferContent);
    return numberOfConsumedBytes;
}

/* receiverFunc     The function that will be called to store the resulting decompressed data */


void inflaterDestroy(Inflater* inf) {
    if (inf) {
        /* delete inf->obj; */
        if (0!=(inf->flags & InfFlags_FreeInputBuffer  )) { free((void*)inf->helperInput.buffer);  }
        if (0!=(inf->flags & InfFlags_FreeOutputBuffer )) { free((void*)inf->helperOutput.buffer); }
        if (0!=(inf->flags & InfFlags_FreeItself       )) { free((void*)inf);                      }
    }
}
