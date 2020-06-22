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


#define inf (*infptr)


//======================================================================================================================
#   pragma mark - READING BIT STREAM

/** Loads the next byte (8 bits) from the input stream to the bitbuffer */
static InfBool InfBS_LoadNextByte(Inflater* infptr) {
    if (inf.inputChunkPtr==inf.inputChunkEnd) { return Inf_FALSE; }
    inf.bitbuffer  |= (*inf.inputChunkPtr++ << inf.bitcounter);
    inf.bitcounter += 8;
    return Inf_TRUE;
}

/** Reads a sequence of bits from the bitbuffer. Returns FALSE if bitbuffer doesn't have enought bits loaded. */
static InfBool InfBS_ReadBits(Inflater* infptr, unsigned* dest, int numberOfBitsToRead) {
    assert( infptr!=NULL && dest!=NULL && numberOfBitsToRead>=0 );
    while (inf.bitcounter<numberOfBitsToRead) { if (!InfBS_LoadNextByte(infptr)) { return Inf_FALSE; }  }
    (*dest) = inf.bitbuffer & ((1<<numberOfBitsToRead)-1);
    inf.bitbuffer  >>= numberOfBitsToRead;
    inf.bitcounter  -= numberOfBitsToRead;
    return Inf_TRUE;
}

/** Reads a sequence of huffman encoded bits from the buffer. Returns FALSE if bitbuffer doesn't have enought bits loaded. */
static InfBool InfBS_ReadCompressedBits(Inflater* infptr, unsigned* dest, const PDZip::ReversedHuffmanDecoder* decoder) {
    assert( infptr!=NULL && dest!=NULL && decoder!=NULL && decoder->isLoaded() );
    InfHuff data;

    // TODO: optimize `readBits` for when there are enough input buffer data
    if ( inf.bitcounter>0 ) {
        data = decoder->decode8(inf.bitbuffer);
        if (data.value.isvalid && data.value.length<=inf.bitcounter) {
            inf.bitbuffer >>= data.value.length;
            inf.bitcounter -= data.value.length;
            (*dest) = data.value.code;
            return Inf_TRUE;
        }
    }
    if ( inf.bitcounter<8 && !InfBS_LoadNextByte(infptr) ) { return Inf_FALSE; }
    data = decoder->decode8(inf.bitbuffer);
    if (data.value.isvalid) {
        inf.bitbuffer  >>= data.value.length;
        inf.bitcounter  -= data.value.length;
        (*dest) = data.value.code;
        return Inf_TRUE;
    }
    if ( inf.bitcounter<16 && !InfBS_LoadNextByte(infptr) ) { return Inf_FALSE; }
    data = decoder->decode16(inf.bitbuffer,data);
    if (data.value.isvalid) {
        inf.bitbuffer >>= data.value.length;
        inf.bitcounter -= data.value.length;
        (*dest) = data.value.code;
        return Inf_TRUE;
    }
    assert( 0 );
    (*dest) = static_cast<unsigned>(-1);
    return Inf_TRUE;
}

/** Reads a WORD (16 bits) from the bitbuffer. Returns FALSE if bitbuffer doesn't have enought bits loaded. */
static InfBool InfBS_ReadWord(Inflater* infptr, unsigned* dest) {
    static const unsigned numberOfBitsToRead = 16;
    const unsigned        bitsToSkip         = (inf.bitcounter%8);
    assert( dest!=NULL && 8*sizeof(*dest)>=numberOfBitsToRead );
    /* skip any padding bits because bitbuffer must be byte aligned before reading a Word */
    inf.bitbuffer  >>= bitsToSkip;
    inf.bitcounter  -= bitsToSkip;
    while ( inf.bitcounter <numberOfBitsToRead ) { if (!InfBS_LoadNextByte(infptr)) { return Inf_FALSE; } }
    assert( inf.bitcounter==numberOfBitsToRead );
    (*dest) = (inf.bitbuffer>>8 & 0x00FF) | (inf.bitbuffer<<8 & 0xFF00);
    inf.bitbuffer = inf.bitcounter = 0;
    return Inf_TRUE;
}

/** Read a DWORD (32 bits) from the bitbuffer. Returns FALSE if bitbuffer doesn't have enought bits loaded. */
static InfBool InfBS_ReadDWord(Inflater* infptr, unsigned* dest) {
    static const unsigned numberOfBitsToRead = 32; // < number of bits in a DWord
    const unsigned        bitsToSkip         = (inf.bitcounter%8);
    assert( dest!=NULL && 8*sizeof(*dest)>=numberOfBitsToRead );
    /* skip any padding bits because bitbuffer must be byte aligned before reading a Word */
    inf.bitbuffer  >>= bitsToSkip;
    inf.bitcounter  -= bitsToSkip;
    while ( inf.bitcounter <numberOfBitsToRead ) { if (!InfBS_LoadNextByte(infptr)) { return Inf_FALSE; } }
    assert( inf.bitcounter==numberOfBitsToRead );
    (*dest) = (inf.bitbuffer>>24 & 0x000000FF) | (inf.bitbuffer>>8 & 0x0000FF00) | (inf.bitbuffer<<8 & 0x00FF0000) | (inf.bitbuffer<<24 & 0xFF000000);
    inf.bitbuffer  = inf.bitcounter = 0;
    return Inf_TRUE;
}

/** Reads a sequence of bytes directly from the input stream (bitbuffer must be empty) */
static InfBool InfBS_ReadBytes(Inflater* infptr, unsigned char* dest, size_t* inout_numberOfBytes) {
    assert( infptr!=NULL && dest!=NULL && inout_numberOfBytes!=NULL && inf.bitcounter==0 );
    const size_t numberOfBytesToRead = (*inout_numberOfBytes);
    const size_t successfulBytes     = min( (inf.inputChunkEnd-inf.inputChunkPtr), numberOfBytesToRead );
    memcpy(dest, inf.inputChunkPtr, successfulBytes);
    inf.inputChunkPtr += successfulBytes;
    (*inout_numberOfBytes) = successfulBytes;
    return (successfulBytes == numberOfBytesToRead);
}


//======================================================================================================================
#   pragma mark - READING CODE/LENGTHS PAIRS

#define InfCL_Add(code, length)                           \
    assert( 0<=code   &&   code<=Inf_LastValidCode   );   \
    assert( 0<=length && length<=Inf_LastValidLength );   \
    assert( inf.cl.insertPtr[length]!=NULL );             \
    if (length>0) {                                       \
        *inf.cl.insertPtr[length]  |= inf.cl.nextIndex;   \
        *( inf.cl.insertPtr[length] = &inf.cl.table[ inf.cl.nextIndex++ ] ) = (code)<<16 | (length)<<10; \
    }                                                     \
    ++inf.cl.size;

#define InfCL_AddArray(array, arraySize) \
    for (int i=0; i<arraySize; i++) {    \
        InfCL_Add(i, array[i]);          \
    }

static void InfCL_AddRange(Inflater* infptr, int firstCode, int lastCode, unsigned length) {
    assert( infptr!=NULL );
    assert( 0<=firstCode && firstCode<=(Inf_LastValidCode)   );
    assert( 0<=lastCode  &&  lastCode<=(Inf_LastValidCode+1) );
    assert( firstCode<=lastCode );
    for (int code=firstCode; code<lastCode; ++code) {
        InfCL_Add(code,length);
    }
}


static void InfCL_Open(Inflater* infptr, InfBool resetRepetitions) {
    if (resetRepetitions) { inf.cl.command = inf.cl.length = inf.cl.repetitions = 0; }
    inf.cl.code      = 0;
    inf.cl.size      = 0;
    inf.cl.nextIndex = (Inf_LastValidLength+1);
    for (int length=0; length<=Inf_LastValidLength; ++length) {
        *( inf.cl.insertPtr[length] = &inf.cl.table[length] ) = 0;
    }
}

static const unsigned* InfCL_Close(Inflater* infptr) {
    int prevLength = 0;
    int length     = 0;
    do {
        unsigned nextIndex = 0;
        do {
            nextIndex = (Inf_NextIndexMask & inf.cl.table[++length]);
        } while ( length<Inf_LastValidLength && nextIndex==0 );
        *inf.cl.insertPtr[prevLength] |= nextIndex;
        prevLength = length;
    } while ( length<Inf_LastValidLength );
    inf.cl.nextIndex = 0;
    return inf.cl.table;
}

static InfBool InfCL_ReadCodes(Inflater* infptr, unsigned numberOfCodes) {
    assert( numberOfCodes>0 );
    assert( infptr!=NULL );
    static const unsigned MaximumNumberOfCodes = 19;
    static const unsigned order[MaximumNumberOfCodes] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };
    
    const unsigned end = (numberOfCodes<MaximumNumberOfCodes) ? numberOfCodes : MaximumNumberOfCodes;
    while (inf.cl.code<end) {
        unsigned length;
        if ( !InfBS_ReadBits(infptr,&length,3) ) { return Inf_FALSE; }
        inf.cl.lengths[ order[inf.cl.code++] ] = length;
    }
    while (inf.cl.code<MaximumNumberOfCodes) { inf.cl.lengths[ order[inf.cl.code++] ] = 0; }
    InfCL_AddArray(inf.cl.lengths, MaximumNumberOfCodes);
    return Inf_TRUE;
}

static InfBool InfCL_ReadCompressedCodes(Inflater* infptr, unsigned numberOfCodes, const PDZip::ReversedHuffmanDecoder* decoder) {
    assert( numberOfCodes>0 );
    assert( infptr!=NULL );
    static const unsigned MaxValidLength = 15; // < the maximum length that can be read
    enum Command {
        Command_ReadNext            = 0,  // < load next command
        Command_CopyPreviousLength  = 16, // < repeat the previous length
        Command_RepeatZeroLength_3  = 17, // < repeat zero length (3 or more times)
        Command_RepeatZeroLength_11 = 18  // < repeat zero length (11 or more times)
    };
    unsigned value = 0;

    // IMPORTANT: add any repetition that is pending from a previous load (ex: literals-table > distance-table)
    while (inf.cl.code<numberOfCodes && inf.cl.repetitions>0) {
        InfCL_Add(inf.cl.code, inf.cl.length);
        ++inf.cl.code; --inf.cl.repetitions;
    }
    
    // add (one by one) all codes with their codeLengths taking into account the number of repetitions
    while ( inf.cl.code<numberOfCodes ) {
        assert( inf.cl.repetitions==0 );
        
        // read a new command
        if ( inf.cl.command==Command_ReadNext ) {
            if ( !InfBS_ReadCompressedBits(infptr, &inf.cl.command, decoder) ) { return Inf_FALSE; }
        }
        // process simple command
        if ( inf.cl.command<=MaxValidLength ) {
            inf.cl.length=inf.cl.command;
            InfCL_Add(inf.cl.code, inf.cl.length);
            ++inf.cl.code;
        }
        // process command with repetition
        else {
            switch (inf.cl.command) {
                case Command_CopyPreviousLength:
                    if (inf.cl.code==0) { /* codeLengths.markAsError(); */ return Inf_TRUE; }
                    if ( !InfBS_ReadBits(infptr,&value,2) ) { return Inf_FALSE; }
                    inf.cl.repetitions = 3 + value;
                    break;
                case Command_RepeatZeroLength_3:
                    if ( !InfBS_ReadBits(infptr,&value,3) ) { return Inf_FALSE; }
                    inf.cl.length      = 0;
                    inf.cl.repetitions = 3 + value;
                    break;
                case Command_RepeatZeroLength_11:
                    if ( !InfBS_ReadBits(infptr,&value,7) ) { return Inf_FALSE; }
                    inf.cl.length      = 0;
                    inf.cl.repetitions = 11 + value;
                    break;
                default:
                    /* codeLengths.markAsError(); */
                    return Inf_TRUE;
            }
            while (inf.cl.code<numberOfCodes && inf.cl.repetitions>0) {
                InfCL_Add(inf.cl.code, inf.cl.length);
                ++inf.cl.code; --inf.cl.repetitions;
            }
        }
        // mark current command as finished
        inf.cl.command = Command_ReadNext;
    }
    return Inf_TRUE;
}



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


/* TODO: remove these globals! */
static PDZip::ReversedHuffmanDecoder  _huffmanDecoder0;
static PDZip::ReversedHuffmanDecoder  _huffmanDecoder1;
static PDZip::ReversedHuffmanDecoder  _huffmanDecoder2;

static PDZip::ReversedHuffmanDecoder* getFixedLiteralDecoder(Inflater* infptr) {
    static PDZip::ReversedHuffmanDecoder decoder;
    if ( !decoder.isLoaded() ) {
        InfCL_Open(infptr, Inf_TRUE);
        InfCL_AddRange(infptr,   0,144, 8);
        InfCL_AddRange(infptr, 144,256, 9);
        InfCL_AddRange(infptr, 256,280, 7);
        InfCL_AddRange(infptr, 280,288, 8);
        decoder.load( InfCL_Close(infptr) );
    }
    return &decoder;
}

static PDZip::ReversedHuffmanDecoder* getFixedDistanceDecoder(Inflater* infptr) {
    static PDZip::ReversedHuffmanDecoder decoder;
    if ( !decoder.isLoaded() ) {
        InfCL_Open(infptr, Inf_TRUE);
        InfCL_AddRange(infptr, 0,32, 5);
        decoder.load( InfCL_Close(infptr) );
    }
    return &decoder;
}



/**
 * Decompress a chunk of compressed data
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
InfAction Inf_decompress(Inflater*                  infptr,
                         unsigned char* const       outputPtr,
                         size_t*                    inout_outputBytes,
                         unsigned char* const       outputBufferBegin,
                         size_t                     outputBufferSize)
{
    assert( inf.inputChunkPtr!=NULL );
    assert( inf.inputChunkPtr<inf.inputChunkEnd );
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

    unsigned char* const writeEnd = outputPtr + (*inout_outputBytes);
    unsigned char*       writePtr = outputPtr;

    unsigned       temp;
    size_t         numberOfBytes;
    bool           canReadAll;
    unsigned char* sequencePtr;
    
    
    InfStep   step        = (InfStep)inf.step;
    InfAction exit_status = (InfAction)(step!=InfStep_End && step!=InfStep_FatalError ? InfAction_ProcessNextChunk : InfError_Failed);
    while ( exit_status==InfAction_ProcessNextChunk ) {
        switch (step) {
                
            case InfStep_Start:
                inf._lastBlock = false;
                op_fallthrough(InfStep_ReadZLibHeader);
                
            //-------------------------------------------------------------------------------------
            // infstep: READ_ZLIB_HEADER
            //
            case InfStep_ReadZLibHeader:
                if ( !InfBS_ReadWord(infptr,&inf._zlibheader) ) { op_return(InfAction_FillInputBuffer); }
                op_fallthrough(InfStep_ReadZLibHeader2);
                
            case InfStep_ReadZLibHeader2:
                if ( (inf._zlibheader % 31)!=0 ) { op_fatal_error(InfError_BadZLibHeader); }
                inf._zlibheader_method = inf._zlibheader>>8 & 0x0F;
                inf._zlibheader_wsize  = 1 << (8 + (inf._zlibheader>>12 & 0x0F));
                inf._zlibheader_level  = inf._zlibheader>>6 & 0x03;
                if (inf._zlibheader_method!=8)       { op_fatal_error(InfError_UnsupportedZLibHeader); }
                if (inf._zlibheader&0x20)            { op_fatal_error(InfError_UnsupportedZLibHeader); }
                if (inf._zlibheader_wsize>(32*1024)) { op_fatal_error(InfError_UnsupportedZLibHeader); }
                //if (_zlibheader_wsize>wrappingBufferSize) { op_fatal_error(Status_UnsupportedZLibHeader); }
                op_fallthrough(InfStep_ProcessNextBlock);
                
            //-------------------------------------------------------------------------------------
            // infstep: READ_BLOCK_HEADER
            //
            case InfStep_ProcessNextBlock:
                if (inf._lastBlock) {
                    printf(" > END OF STREAM\n\n");
                    op_end(InfAction_Finish);
                }
                op_fallthrough(InfStep_ReadBlockHeader);
                
            case InfStep_ReadBlockHeader:
                if ( !InfBS_ReadBits(infptr,&temp,3) ) { op_return(InfAction_FillInputBuffer); }
                inf._lastBlock = (temp&0x01)==0x01;
                inf._blocktype = (temp>>1);
                if      (inf._blocktype==BlockType_Uncompressed)   {
                    printf(" > Start Uncompressed Block\n");
                    op_goto(InfStep_Process_UncompressedBlock);
                }
                else if (inf._blocktype==BlockType_FixedHuffman)   {
                    printf(" > Start Fixed Huffman Block\n");
                    op_goto(InfStep_Load_FixedHuffmanDecoders);
                }
                else if (inf._blocktype==BlockType_DynamicHuffman) {
                    printf(" > Start Dynamic Huffman Block\n");
                    op_goto(InfStep_Load_DynamicHuffmanDecoders);
                }
                printf(" > Fatal Error (Bad Block Type)\n");
                op_fatal_error(InfError_BadBlockType);
                
            //-------------------------------------------------------------------------------------
            // Process_UncompressedBlock
            //
            case InfStep_Process_UncompressedBlock:
                if ( !InfBS_ReadDWord(infptr,&temp) ) { op_return(InfAction_FillInputBuffer); }
                inf._seq_len = (temp >> 16);
                if ( inf._seq_len != ((~temp)&0xFFFF) ) { op_fatal_error(InfError_BadBlockLength); }
                op_fallthrough(InfStep_Output_UncompressedBlock);
                
            case InfStep_Output_UncompressedBlock:
                numberOfBytes = min( (writeEnd-writePtr), inf._seq_len );
                canReadAll    = InfBS_ReadBytes(infptr, writePtr, &numberOfBytes);
                inf._seq_len -= numberOfBytes;
                writePtr     += numberOfBytes;
                if ( inf._seq_len>0 )   { op_return(InfAction_UseOutputBufferContent); }
                else if ( !canReadAll ) { op_return(InfAction_FillInputBuffer);        }
                
                assert( inf._seq_len==0 );
                op_goto(InfStep_ProcessNextBlock);

            //-------------------------------------------------------------------------------------
            // Load FixedHuffmanDecoders
            //
            case InfStep_Load_FixedHuffmanDecoders:
                inf._literalDecoder  = getFixedLiteralDecoder(infptr);
                inf._distanceDecoder = getFixedDistanceDecoder(infptr);
                op_goto(InfStep_Process_CompressedBlock);
                
            //-------------------------------------------------------------------------------------
            // Load DynamicHuffmanDecoders
            //
            case InfStep_Load_DynamicHuffmanDecoders:
                if ( !InfBS_ReadBits(infptr,&temp,5+5+4) ) { op_return(InfAction_FillInputBuffer); }
                inf._hlit  = temp     & 0x1F;
                inf._hdist = temp>>5  & 0x1F;
                inf._hclen = temp>>10 & 0x0F;
                op_fallthrough(InfStep_Load_FrontHuffmanTable);

            //----------------------------------
            // load "front" huffman table
            //
            case InfStep_Load_FrontHuffmanTable:
                InfCL_Open(infptr,Inf_TRUE);
                op_fallthrough(InfStep_Load_FrontHuffmanTable2);
                
            case InfStep_Load_FrontHuffmanTable2:
                if ( !InfCL_ReadCodes(infptr,inf._hclen+4) ) { op_return(InfAction_FillInputBuffer); }
                inf._frontDecoder = &_huffmanDecoder0;
                inf._frontDecoder->load( InfCL_Close(infptr) );
                op_fallthrough(InfStep_Load_LiteralHuffmanTable);
                
            //----------------------------------
            // load literal-length huffman table
            //
            case InfStep_Load_LiteralHuffmanTable:
                InfCL_Open(infptr,Inf_TRUE);
                op_fallthrough(InfStep_Load_LiteralHuffmanTable2);
                
            case InfStep_Load_LiteralHuffmanTable2:
                if ( !InfCL_ReadCompressedCodes(infptr,inf._hlit+257,inf._frontDecoder) ) { op_return(InfAction_FillInputBuffer); }
                inf._literalDecoder = &_huffmanDecoder1;
                inf._literalDecoder->load( InfCL_Close(infptr) );
                op_fallthrough(InfStep_Load_DistanceHuffmanTable);
                
            //----------------------------------
            // load distance huffman table
            //
            case InfStep_Load_DistanceHuffmanTable:
                InfCL_Open(infptr,Inf_FALSE);
                op_fallthrough(InfStep_Load_DistanceHuffmanTable2);
                
            case InfStep_Load_DistanceHuffmanTable2:
                if ( !InfCL_ReadCompressedCodes(infptr,inf._hdist+1,inf._frontDecoder) ) { op_return(InfAction_FillInputBuffer); }
                inf._distanceDecoder = &_huffmanDecoder2;
                inf._distanceDecoder->load( InfCL_Close(infptr) );
                op_fallthrough(InfStep_Process_CompressedBlock);
                
            //-------------------------------------------------------------------------------------
            // Process Compressed Block
            //
            case InfStep_Process_CompressedBlock:
            case InfStep_Read_LiteralOrLength:
                if ( !InfBS_ReadCompressedBits(infptr,&inf._literal,inf._literalDecoder) ) { op_return(InfAction_FillInputBuffer); }
                op_fallthrough(InfStep_Read_LiteralOrLength2);
                
            case InfStep_Read_LiteralOrLength2:
                if (inf._literal <Inf_EndOfBlock) {
                    if (writePtr==writeEnd) { op_return(InfAction_UseOutputBufferContent); }
                    *writePtr++ = inf._literal;
                    op_goto(InfStep_Read_LiteralOrLength);
                }
                else if (inf._literal==Inf_EndOfBlock) {
                    printf(" > EndOfBlock\n");
                    op_goto(InfStep_ProcessNextBlock);
                }
                else if (inf._literal>Inf_MaxValidLengthCode) { op_fatal_error(InfError_BadBlockContent); }
                inf._literal -= 257;
                op_fallthrough(InfStep_Read_LengthBits);
                
            case InfStep_Read_LengthBits:
                if ( !InfBS_ReadBits(infptr, &temp, lengthExtraBits[inf._literal]) ) { op_return(InfAction_FillInputBuffer); }
                inf._seq_len   =   temp + lengthStarts[inf._literal];
                assert( inf._seq_len<(32*1024) );
                op_fallthrough(InfStep_Read_Distance);
                
            case InfStep_Read_Distance:
                if ( !InfBS_ReadCompressedBits(infptr,&inf._literal,inf._distanceDecoder) ) { op_return(InfAction_FillInputBuffer); }
                if (inf._literal>Inf_MaxValidDistanceCode) { op_fatal_error(InfError_BadBlockContent); }
                op_fallthrough(InfStep_Read_DistanceBits);
                
            case InfStep_Read_DistanceBits:
                if ( !InfBS_ReadBits(infptr, &temp, distanceExtraBits[inf._literal]) ) { op_return(InfAction_FillInputBuffer); }
                inf._seq_dist =  temp + distanceStarts[inf._literal];
                assert( inf._seq_dist<(32*1024) );
                op_fallthrough(InfStep_Output_RepeatedSequence);

            //-------------------------------------------------------------------------------------
            // Output_RepeatedSequence
            //
            //         outputBegin     writePtr                      writeEnd
            //              V             V                             V
            // <... ghost > # [[ straight * overlapped  ...... ghost ]] #
            //
            case InfStep_Output_RepeatedSequence:
                sequencePtr = writePtr - inf._seq_dist;
                if (sequencePtr<outputBufferBegin) {
                    sequencePtr += outputBufferSize;
                    //-- copy bytes ---
                    inf._seq_len -= numberOfBytes = min( inf._seq_len, (writeEnd-sequencePtr) );
                    memmove( writePtr, sequencePtr, numberOfBytes ); writePtr+=numberOfBytes;
                    //-----------------
                    if ( inf._seq_len==0 ) { op_goto(InfStep_Read_LiteralOrLength); }
                    if ( writePtr==writeEnd )  { op_return(InfAction_UseOutputBufferContent); }
                    sequencePtr = outputBufferBegin;
                }
                //-- copy bytes ---
                inf._seq_len -= numberOfBytes = min( inf._seq_len, (writeEnd-writePtr) );
                while ( numberOfBytes-->0 ) { *writePtr++ = *sequencePtr++; }
                //-----------------
                if ( inf._seq_len==0 ) { op_goto(InfStep_Read_LiteralOrLength); }
                if ( writePtr==writeEnd )  { op_return(InfAction_UseOutputBufferContent); }
                op_fatal_error(InfError_BadBlockContent);
                
                
            case InfStep_FatalError:
            case InfStep_End:
                break;
        }
    }
    
    
    inf.step                     = step;
    inf.outputChunkSize          = (writePtr          - outputPtr);
    inf.outputBufferContentSize += inf.outputChunkSize;
    inf.action                   = exit_status>=0 ? (InfAction)exit_status : InfAction_Finish;
    inf.error                    = exit_status <0 ? (InfError )exit_status : InfError_None;
    
    printf(" # exit_stats = %d\n", exit_status);
    return exit_status;
}

#undef inf




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
    
    inf->inputChunkPtr   = NULL;
    inf->inputChunkEnd   = NULL;
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
    inf->step       = 0;
    inf->bitbuffer  = 0;
    inf->bitcounter = 0;
    
    
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
    size_t outputBytes;
    /* const Byte* const inputBufferEnd  = (Byte*)inputBuffer  + inputBufferSize; */
    Byte* const       outputBufferEnd = (Byte*)outputBuffer + outputBufferSize;
    
    assert( inputBuffer !=NULL && inputBufferSize >0 );
    assert( outputBuffer!=NULL && outputBufferSize>0 );
    
    if (inflater->finished) {
        inflater->action=InfAction_Finish;
    }
    
    switch ( inflater->action ) {
            
        case InfAction_Init:
            inflater->inputChunkPtr = (const Byte*)inputBuffer;
            inflater->inputChunkEnd = inflater->inputChunkPtr + inputBufferSize;
            inflater->outputChunk   = (Byte      *)outputBuffer;
            break;
            
        case InfAction_FillInputBuffer:
            inflater->inputChunkPtr = (const Byte*)inputBuffer;
            inflater->inputChunkEnd = inflater->inputChunkPtr + inputBufferSize;
            inflater->outputChunk += inflater->outputChunkSize;
            break;
            
        case InfAction_UseOutputBufferContent:
            inflater->outputChunk  = (Byte*)outputBuffer;
            inflater->outputBufferContentSize = 0;
            break;
            
        case InfAction_ProcessNextChunk:
            inflater->outputChunk += inflater->outputChunkSize;
            break;
            
        /* InfAction_Finish */
        default:
            inflater->outputChunkSize = 0;
            return inflater->action;
    }
    
    assert( inflater->inputChunkPtr < inflater->inputChunkEnd  );
    assert( inflater->outputChunk   < outputBufferEnd          );
    
    /* inputBytes  = (inputBufferEnd  - inflater->inputChunk ); */
    outputBytes = (outputBufferEnd - inflater->outputChunk);
    
    status = Inf_decompress(inflater,
                            inflater->outputChunk, &outputBytes,
                            (unsigned char*)outputBuffer, outputBufferSize);
    
    
    /* hack to add an 'UseOutputBufferContent' before 'Finish' */
    inflater->finished = (inflater->action==InfAction_Finish);
    if (inflater->finished && inflater->outputChunkSize>0) {
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
    InfAction action; /* size_t numberOfConsumedBytes=0; */
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
        /* numberOfConsumedBytes += inf->inputChunkSize; */
    } while (action==InfAction_ProcessNextChunk || action==InfAction_UseOutputBufferContent);
    
/*  return numberOfConsumedBytes;  */
    return compressedDataSize;
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
