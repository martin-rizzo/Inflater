/**
 * @file       inflater.c
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
#include <assert.h>
#include <stdlib.h>
#include <string.h>
/* #include "PDZip_ReversedHuffmanDecoder.h" */

#include <stdio.h>

static size_t min(size_t a, size_t b) { return a<b ? a : b; }


#define Inf_MainTableSize 256      /**< 8bits                        */
#define InfHuffTableSize  (2*1024) /**< main-table + all sub-tables  */
#define Inf_LastValidLength 18
#define Inf_InvalidLength   24
#define Inf_MaxNumberOfCodes 19
#define Inf_CodeLengthTableSize ((Inf_LastValidLength+1)+(Inf_LastValidCode+1))


#define inf (*infptr)


/*====================================================================================================================*/
#   pragma mark - HUFFMAN DECODER

static int CL_getFirstIndex(const unsigned* sourtable) {
    assert( sourtable!=NULL );
    return 0x03FF & sourtable[0];
}

static int CL_getNextIndex(const unsigned* sourtable, int index) {
    assert( sourtable!=NULL );
    assert( 0<=index && index<Inf_CodeLengthTableSize );
    return 0x03FF & sourtable[index];
}


static unsigned CL_getCodeAt(const unsigned* sourtable, int index) {
    assert( sourtable!=NULL );
    assert( 0<=index && index<Inf_CodeLengthTableSize );
    return 0xFFFF & sourtable[index]>>16;
}

static unsigned CL_getLengthAt(const unsigned* sourtable, int index) {
    assert( sourtable!=NULL );
    assert( 0<=index && index<Inf_CodeLengthTableSize );
    return 0x3F & sourtable[index]>>10;
}

static unsigned reversedINC(unsigned huffman, unsigned length) {
    unsigned incBit;
    assert( length>0 );
    incBit = 1<<(length-1); do { huffman ^= incBit; } while ( (huffman&incBit)==0 && (incBit>>=1)!=0 );
    return huffman;
}

static void setWithRepetitions(InfHuff*  table,
                               unsigned  availableEntries,
                               unsigned  huffman,
                               InfHuff   data,
                               unsigned  length)
{
    unsigned unknownBits; const unsigned unknownStep = (1<<length);
    assert( table!=NULL && huffman<(1<<length) && length>=1 );
    for (unknownBits=0; unknownBits<availableEntries; unknownBits+=unknownStep) {
        table[ unknownBits | huffman ] = data;
    }
}

static int CL_makeSubTable(InfHuff*           subtable,
                           int                availableEntries,
                           unsigned           firstHuffman,
                           unsigned           maxLength,
                           const unsigned*    sourtable,
                           int                firstIndex,
                           int                endIndex)
{
    static const unsigned DiscardedBits = 8; /**< number of bits discarded by the subtable */
    const int numberOfEntries = 1<<(maxLength-DiscardedBits);
    unsigned huffman; int index;
    
    assert( numberOfEntries<=availableEntries  );
    
    huffman = firstHuffman;
    index   = firstIndex;
    while ( index!=endIndex ) {
        const unsigned code   = CL_getCodeAt(sourtable,index);
        const unsigned length = CL_getLengthAt(sourtable,index);
        InfHuff data;
        data.value.code    = code;
        data.value.length  = length;
        data.value.isvalid = 1;
        setWithRepetitions(subtable, numberOfEntries, huffman>>DiscardedBits, data, length-DiscardedBits);

        huffman = reversedINC(huffman,length);
        index   = CL_getNextIndex(sourtable,index);
    }
    return numberOfEntries;
}


const InfHuff* InfHD_load(InfHuff* table, const unsigned int *sourtable) {
    InfHuff invalid, data; int i, index, insertIndex;
    unsigned huffman, length, maxLength;
    assert( table!=NULL && sourtable!=NULL  );
    
    /* reset the main-table */
    invalid.value.isvalid = 1;
    invalid.value.code    = 0;
    invalid.value.length  = Inf_InvalidLength;
    for (i=0; i<Inf_MainTableSize; ++i) { table[i] = invalid; }
    
    huffman = 0;
    index   = CL_getFirstIndex(sourtable);
    length  = CL_getLengthAt(sourtable,index);

    /* lengths from 1 to 8                              */
    /* unknown bits are filled with all possible values */
    while ( index!=0 && length<=8 ) {
        InfHuff data;
        data.value.isvalid = 1;
        data.value.code    = CL_getCodeAt(sourtable,index);
        data.value.length  = length;
        setWithRepetitions(table, 256, huffman, data, length);
        huffman = reversedINC(huffman,length);
        index   = CL_getNextIndex(sourtable,index);
        length  = CL_getLengthAt(sourtable,index);
    }
    /* lengths from 9         */
    /* subtables are created  */
    insertIndex = Inf_MainTableSize;
    while ( index!=0 ) {
        const unsigned firstHuffman = huffman;
        const int      firstIndex   = index;
        do {
            length  = CL_getLengthAt(sourtable,index);
            huffman = reversedINC(huffman,length);
            index   = CL_getNextIndex(sourtable,index);
        } while ( index!=0 && (huffman&0xFF)==(firstHuffman&0xFF) );
        
        maxLength = length;
        data.subtable.error = 0;
        data.subtable.index = insertIndex;
        data.subtable.mask  = (1<<(maxLength-8))-1;
        table[firstHuffman] = data;
        insertIndex += CL_makeSubTable(&table[insertIndex],
                                       (InfHuffTableSize-insertIndex),
                                       firstHuffman, maxLength,
                                       sourtable, firstIndex, index
                                       );
    }
    return table;
}


#define InfHD_decode(table, tmp, bitbuffer) \
    (tmp=table[bitbuffer & 0xFF], tmp.value.isvalid ? tmp : table[ tmp.subtable.index + ((bitbuffer>>8) & tmp.subtable.mask) ])


/*====================================================================================================================*/
#   pragma mark - BIT STREAM

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
    
    while (inf.bitcounter<numberOfBitsToRead) {
        if (!InfBS_LoadNextByte(infptr)) { return Inf_FALSE; }
    }
    (*dest) = inf.bitbuffer & ((1<<numberOfBitsToRead)-1);
    inf.bitbuffer  >>= numberOfBitsToRead;
    inf.bitcounter  -= numberOfBitsToRead;
    return Inf_TRUE;
}

/** Reads a sequence of huffman encoded bits from the buffer. Returns FALSE if bitbuffer doesn't have enought bits loaded. */
static InfBool InfBS_ReadCompressedBits(Inflater* infptr, unsigned* dest, const InfHuff* table) {
    unsigned numberOfBitsToAdvance; InfHuff data, tmp;
    assert( infptr!=NULL && dest!=NULL && table!=NULL );

    if (inf.bitcounter==0 && !InfBS_LoadNextByte(infptr)) { return Inf_FALSE; }
    data=InfHD_decode(table, tmp, inf.bitbuffer); numberOfBitsToAdvance=data.value.length;
    
    while (inf.bitcounter<numberOfBitsToAdvance) {
        if (!InfBS_LoadNextByte(infptr)) { return Inf_FALSE; }
        data=InfHD_decode(table, tmp, inf.bitbuffer); numberOfBitsToAdvance=data.value.length;
    }
    (*dest) = data.value.code;
    inf.bitbuffer >>= numberOfBitsToAdvance;
    inf.bitcounter -= numberOfBitsToAdvance;
    return Inf_TRUE;
}

/** Read a DWORD (32 bits) from the bitbuffer. Returns FALSE if bitbuffer doesn't have enought bits loaded. */
static InfBool InfBS_ReadDWord(Inflater* infptr, unsigned* dest) {
    static const unsigned numberOfBitsToRead = 32; /**< number of bits in a DWord */
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
    size_t numberOfBytesToRead, successfulBytes;
    assert( infptr!=NULL && dest!=NULL && inout_numberOfBytes!=NULL && inf.bitcounter==0 );
    numberOfBytesToRead = (*inout_numberOfBytes);
    successfulBytes     = min( (inf.inputChunkEnd-inf.inputChunkPtr), numberOfBytesToRead );
    memcpy(dest, inf.inputChunkPtr, successfulBytes);
    inf.inputChunkPtr += successfulBytes;
    (*inout_numberOfBytes) = successfulBytes;
    return (successfulBytes == numberOfBytesToRead);
}


/*====================================================================================================================*/
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

static void InfCL_AddRange(Inflater* infptr, int firstCode, int lastCode, unsigned length) {
    int code;
    assert( infptr!=NULL );
    assert( 0<=firstCode && firstCode<=(Inf_LastValidCode)   );
    assert( 0<=lastCode  &&  lastCode<=(Inf_LastValidCode+1) );
    assert( firstCode<=lastCode );
    for (code=firstCode; code<lastCode; ++code) {
        InfCL_Add(code,length);
    }
}


static void InfCL_Open(Inflater* infptr, InfBool resetRepetitions) {
    int length;
    if (resetRepetitions) { inf.cl.command = inf.cl.length = inf.cl.repetitions = 0; }
    inf.cl.code      = 0;
    inf.cl.size      = 0;
    inf.cl.nextIndex = (Inf_LastValidLength+1);
    for (length=0; length<=Inf_LastValidLength; ++length) {
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
    static const unsigned order[Inf_MaxNumberOfCodes] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };
    const unsigned end = (numberOfCodes<Inf_MaxNumberOfCodes) ? numberOfCodes : Inf_MaxNumberOfCodes;
    int i;
    assert( infptr!=NULL && numberOfCodes>0 );
    
    while (inf.cl.code<end) {
        unsigned length;
        if ( !InfBS_ReadBits(infptr,&length,3) ) { return Inf_FALSE; }
        inf.cl.lengths[ order[inf.cl.code++] ] = length;
    }
    while (inf.cl.code<Inf_MaxNumberOfCodes) { inf.cl.lengths[ order[inf.cl.code++] ] = 0; }
    for (i=0; i<Inf_MaxNumberOfCodes; i++) { InfCL_Add(i, inf.cl.lengths[i]);  }
    return Inf_TRUE;
}

static InfBool InfCL_ReadCompressedCodes(Inflater* infptr, unsigned numberOfCodes,  const InfHuff* table) {
    enum Command {
        Command_ReadNext            = 0,  /**< load next command                             */
        Command_UseCommandAsLength  = 15, /**< length is the value stored in 'infcl.command' */
        Command_CopyPreviousLength  = 16, /**< repeat the previous length                    */
        Command_RepeatZeroLength_3  = 17, /**< repeat zero length (3 or more times)          */
        Command_RepeatZeroLength_11 = 18  /**< repeat zero length (11 or more times)         */
    };
    unsigned value = 0;
    assert( infptr!=NULL && numberOfCodes>0 && table!=NULL );

    /* IMPORTANT: add any repetition that is pending from a previous load (ex: literals-table > distance-table) */
    while (inf.cl.code<numberOfCodes && inf.cl.repetitions>0) {
        InfCL_Add(inf.cl.code, inf.cl.length);
        ++inf.cl.code; --inf.cl.repetitions;
    }
    
    /* add (one by one) all codes with their codeLengths taking into account the number of repetitions */
    while ( inf.cl.code<numberOfCodes ) {
        assert( inf.cl.repetitions==0 );
        
        /* read a new command */
        if ( inf.cl.command==Command_ReadNext ) {
            if ( !InfBS_ReadCompressedBits(infptr, &inf.cl.command, table) ) { return Inf_FALSE; }
        }
        /* process simple command */
        if ( inf.cl.command<=Command_UseCommandAsLength ) {
            inf.cl.length=inf.cl.command;
            InfCL_Add(inf.cl.code, inf.cl.length);
            ++inf.cl.code;
        }
        /* process command with repetition */
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
        /* mark current command as finished */
        inf.cl.command = Command_ReadNext;
    }
    return Inf_TRUE;
}



/*====================================================================================================================*/
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
    BlockType_Uncompressed   = 0,
    BlockType_FixedHuffman   = 1,
    BlockType_DynamicHuffman = 2
} BlockType;




typedef enum InfStep {
    InfStep_Start,
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
static InfHuff huffmanTable0[InfHuffTableSize];
static InfHuff huffmanTable1[InfHuffTableSize];
static InfHuff huffmanTable2[InfHuffTableSize];


static InfHuff* getFixedLiteralDecoder(Inflater* infptr) {
    static InfHuff table[InfHuffTableSize];
    static InfBool loaded = Inf_FALSE;
    if (!loaded) {
        loaded = Inf_TRUE;
        InfCL_Open(infptr, Inf_TRUE);
        InfCL_AddRange(infptr,   0,144, 8);
        InfCL_AddRange(infptr, 144,256, 9);
        InfCL_AddRange(infptr, 256,280, 7);
        InfCL_AddRange(infptr, 280,288, 8);
        InfHD_load(table, InfCL_Close(infptr) );
    }
    return table;
}

static InfHuff* getFixedDistanceDecoder(Inflater* infptr) {
    static InfHuff table[InfHuffTableSize];
    static InfBool loaded = Inf_FALSE;
    if (!loaded) {
        loaded = Inf_TRUE;
        InfCL_Open(infptr, Inf_TRUE);
        InfCL_AddRange(infptr, 0,32, 5);
        InfHD_load(table, InfCL_Close(infptr) );
    }
    return table;
}


#define inf__FILL_INPUT_BUFFER()                             inf.action=InfAction_FillInputBuffer;        break;
#define inf__USE_OUTPUT_BUFFER_CONTENT()                     inf.action=InfAction_UseOutputBufferContent; break;
#define inf__FINISH()               step=InfStep_End;        inf.action=InfAction_Finish;                 break;
#define inf__ERROR(err)             step=InfStep_FatalError; inf.action=InfAction_Finish; inf.error=err;  break;
#define inf__goto(dest_step)        step=dest_step;                                                       break;
#define inf__fallthrough(next_step) step=next_step;


/**
 * Decompress a chunk of compressed data
 *
 * @param outputBufferBegin
 *     The beginning of the output buffer (used when copying previous sequences)
 *
 * @param outputBufferSize
 *     The total size of the output buffer (in number of bytes)
 */
InfAction Inf_decompress(Inflater*                  infptr,
                         unsigned char* const       outputBufferBegin,
                         size_t                     outputBufferSize)
{
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

    unsigned       temp;
    size_t         numberOfBytes;
    InfBool        canReadAll;
    unsigned char* sequencePtr;
    InfStep step;

    unsigned char* const writeEnd = (outputBufferBegin + outputBufferSize);
    unsigned char*       writePtr = inf.outputChunk;

    assert( inf.inputChunkPtr!=NULL );
    assert( inf.inputChunkPtr<inf.inputChunkEnd );
    assert( inf.outputChunk!=NULL && inf.outputChunk<(outputBufferBegin+outputBufferSize) );
    assert( outputBufferBegin!=NULL );
    assert( outputBufferSize>=(32*1024) );
    
    step = (InfStep)inf.step;
    if (step==InfStep_End || step==InfStep_FatalError) { return inf.action; }
    
    inf.action=InfAction_ProcessNextChunk;
    while ( inf.action==InfAction_ProcessNextChunk ) {
        switch (step) {
                
            case InfStep_Start:
                inf._lastBlock = Inf_FALSE;
                inf__fallthrough(InfStep_ProcessNextBlock);
                
            /*-------------------------------------------------------------------------------------
             * infstep: READ_BLOCK_HEADER
             */
            case InfStep_ProcessNextBlock:
                if (inf._lastBlock) {
                    if (writePtr > inf.outputChunk) {
                        printf(" > END (use remaining output buffer)\n");
                        inf__USE_OUTPUT_BUFFER_CONTENT();
                    }
                    printf(" > END OF STREAM\n\n");
                    inf__FINISH();
                }
                inf__fallthrough(InfStep_ReadBlockHeader);
                
            case InfStep_ReadBlockHeader:
                if ( !InfBS_ReadBits(infptr,&temp,3) ) { inf__FILL_INPUT_BUFFER(); }
                inf._lastBlock = (temp&0x01)==0x01;
                inf._blocktype = (temp>>1);
                if      (inf._blocktype==BlockType_Uncompressed)   {
                    printf(" > Start Uncompressed Block\n");
                    inf__goto(InfStep_Process_UncompressedBlock);
                }
                else if (inf._blocktype==BlockType_FixedHuffman)   {
                    printf(" > Start Fixed Huffman Block\n");
                    inf__goto(InfStep_Load_FixedHuffmanDecoders);
                }
                else if (inf._blocktype==BlockType_DynamicHuffman) {
                    printf(" > Start Dynamic Huffman Block\n");
                    inf__goto(InfStep_Load_DynamicHuffmanDecoders);
                }
                printf(" > Fatal Error (Bad Block Type)\n");
                inf__ERROR(InfError_BadBlockType);
                
            /*-------------------------------------------------------------------------------------
             * Process_UncompressedBlock
             */
            case InfStep_Process_UncompressedBlock:
                if ( !InfBS_ReadDWord(infptr,&temp) ) { inf__FILL_INPUT_BUFFER(); }
                inf._seq_len = (temp >> 16);
                if ( inf._seq_len != ((~temp)&0xFFFF) ) { inf__ERROR(InfError_BadBlockLength); }
                inf__fallthrough(InfStep_Output_UncompressedBlock);
                
            case InfStep_Output_UncompressedBlock:
                numberOfBytes = min( (writeEnd-writePtr), inf._seq_len );
                canReadAll    = InfBS_ReadBytes(infptr, writePtr, &numberOfBytes);
                inf._seq_len -= numberOfBytes;
                writePtr     += numberOfBytes;
                if ( inf._seq_len>0 )   { inf__USE_OUTPUT_BUFFER_CONTENT(); }
                else if ( !canReadAll ) { inf__FILL_INPUT_BUFFER();         }
                
                assert( inf._seq_len==0 );
                inf__goto(InfStep_ProcessNextBlock);

            /*-------------------------------------------------------------------------------------
             * Load FixedHuffmanDecoders
             */
            case InfStep_Load_FixedHuffmanDecoders:
                inf.literalDecoder  = getFixedLiteralDecoder (infptr);
                inf.distanceDecoder = getFixedDistanceDecoder(infptr);
                inf__goto(InfStep_Process_CompressedBlock);
                
            /*-------------------------------------------------------------------------------------
             * Load DynamicHuffmanDecoders
             */
            case InfStep_Load_DynamicHuffmanDecoders:
                if ( !InfBS_ReadBits(infptr,&temp,5+5+4) ) { inf__FILL_INPUT_BUFFER(); }
                inf._hlit  = temp     & 0x1F;
                inf._hdist = temp>>5  & 0x1F;
                inf._hclen = temp>>10 & 0x0F;
                inf__fallthrough(InfStep_Load_FrontHuffmanTable);

            /*----------------------------------
             * load "front" huffman table
             */
            case InfStep_Load_FrontHuffmanTable:
                InfCL_Open(infptr,Inf_TRUE);
                inf__fallthrough(InfStep_Load_FrontHuffmanTable2);
                
            case InfStep_Load_FrontHuffmanTable2:
                if ( !InfCL_ReadCodes(infptr,inf._hclen+4) ) { inf__FILL_INPUT_BUFFER(); }
                inf.frontDecoder = InfHD_load(huffmanTable0, InfCL_Close(infptr));
                inf__fallthrough(InfStep_Load_LiteralHuffmanTable);
                
            /*----------------------------------
             * load literal-length huffman table
             */
            case InfStep_Load_LiteralHuffmanTable:
                InfCL_Open(infptr,Inf_TRUE);
                inf__fallthrough(InfStep_Load_LiteralHuffmanTable2);
                
            case InfStep_Load_LiteralHuffmanTable2:
                if ( !InfCL_ReadCompressedCodes(infptr,inf._hlit+257,inf.frontDecoder) ) { inf__FILL_INPUT_BUFFER(); }
                inf.literalDecoder = InfHD_load(huffmanTable1, InfCL_Close(infptr));
                inf__fallthrough(InfStep_Load_DistanceHuffmanTable);
                
            /*----------------------------------
             * load distance huffman table
             */
            case InfStep_Load_DistanceHuffmanTable:
                InfCL_Open(infptr,Inf_FALSE);
                inf__fallthrough(InfStep_Load_DistanceHuffmanTable2);
                
            case InfStep_Load_DistanceHuffmanTable2:
                if ( !InfCL_ReadCompressedCodes(infptr,inf._hdist+1,inf.frontDecoder) ) { inf__FILL_INPUT_BUFFER(); }
                inf.distanceDecoder = InfHD_load(huffmanTable2, InfCL_Close(infptr));
                inf__fallthrough(InfStep_Process_CompressedBlock);
                
            /*-------------------------------------------------------------------------------------
             * Process Compressed Block
             */
            case InfStep_Process_CompressedBlock:
            case InfStep_Read_LiteralOrLength:
                if ( !InfBS_ReadCompressedBits(infptr,&inf._literal,inf.literalDecoder) ) { inf__FILL_INPUT_BUFFER(); }
                inf__fallthrough(InfStep_Read_LiteralOrLength2);
                
            case InfStep_Read_LiteralOrLength2:
                if (inf._literal <Inf_EndOfBlock) {
                    if (writePtr==writeEnd) { inf__USE_OUTPUT_BUFFER_CONTENT(); }
                    *writePtr++ = inf._literal;
                    inf__goto(InfStep_Read_LiteralOrLength);
                }
                else if (inf._literal==Inf_EndOfBlock) {
                    printf(" > EndOfBlock\n");
                    inf__goto(InfStep_ProcessNextBlock);
                }
                else if (inf._literal>Inf_MaxValidLengthCode) { inf__ERROR(InfError_BadBlockContent); }
                inf._literal -= 257;
                inf__fallthrough(InfStep_Read_LengthBits);
                
            case InfStep_Read_LengthBits:
                if ( !InfBS_ReadBits(infptr, &temp, lengthExtraBits[inf._literal]) ) { inf__FILL_INPUT_BUFFER(); }
                inf._seq_len   =   temp + lengthStarts[inf._literal];
                assert( inf._seq_len<(32*1024) );
                inf__fallthrough(InfStep_Read_Distance);
                
            case InfStep_Read_Distance:
                if ( !InfBS_ReadCompressedBits(infptr,&inf._literal,inf.distanceDecoder) ) { inf__FILL_INPUT_BUFFER(); }
                if (inf._literal>Inf_MaxValidDistanceCode) { inf__ERROR(InfError_BadBlockContent); }
                inf__fallthrough(InfStep_Read_DistanceBits);
                
            case InfStep_Read_DistanceBits:
                if ( !InfBS_ReadBits(infptr, &temp, distanceExtraBits[inf._literal]) ) { inf__FILL_INPUT_BUFFER(); }
                inf._seq_dist =  temp + distanceStarts[inf._literal];
                assert( inf._seq_dist<(32*1024) );
                inf__fallthrough(InfStep_Output_RepeatedSequence);

            /*-------------------------------------------------------------------------------------
             * Output_RepeatedSequence
             *
             *         outputBegin     writePtr                      writeEnd
             *              V             V                             V
             * <... ghost > # [[ straight * overlapped  ...... ghost ]] #
             */
            case InfStep_Output_RepeatedSequence:
                sequencePtr = writePtr - inf._seq_dist;
                if (sequencePtr<outputBufferBegin) {
                    sequencePtr += outputBufferSize;
                    /*-- copy bytes ---*/
                    inf._seq_len -= numberOfBytes = min( inf._seq_len, (writeEnd-sequencePtr) );
                    memmove( writePtr, sequencePtr, numberOfBytes ); writePtr+=numberOfBytes;
                    /*-----------------*/
                    if ( inf._seq_len==0    ) { inf__goto(InfStep_Read_LiteralOrLength); }
                    if ( writePtr==writeEnd ) { inf__USE_OUTPUT_BUFFER_CONTENT();        }
                    sequencePtr = outputBufferBegin;
                }
                /*-- copy bytes ---*/
                inf._seq_len -= numberOfBytes = min( inf._seq_len, (writeEnd-writePtr) );
                while ( numberOfBytes-->0 ) { *writePtr++ = *sequencePtr++; }
                /*-----------------*/
                if ( inf._seq_len==0    ) { inf__goto(InfStep_Read_LiteralOrLength); }
                if ( writePtr==writeEnd ) { inf__USE_OUTPUT_BUFFER_CONTENT();        }
                inf__ERROR(InfError_BadBlockContent);
                
                
            case InfStep_FatalError:
            case InfStep_End:
                break;
        }
    }
    
    inf.step                     = (unsigned)step;
    inf.outputChunkSize          = (writePtr - inf.outputChunk);
    inf.outputBufferContentSize += inf.outputChunkSize;
    printf(" # inf.action = %d\n", inf.action);
    return inf.action;
}

#undef inf__FILL_INPUT_BUFFER
#undef inf__USE_OUTPUT_BUFFER_CONTENT
#undef inf__FINISH
#undef inf__ERROR
#undef inf__goto
#undef inf__fallthrough

#undef inf




Inflater* inflaterCreate(void* workingBuffer, size_t workingBufferSize) {
    Inflater* inf = (Inflater*)malloc(sizeof(Inflater));
    
    inf->mode            = InfMode_Uninitialized;
    inf->flags           = InfFlags_FreeItself;
    inf->action          = InfAction_Init;
    inf->error           = InfError_None;
    
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

    assert( inputBuffer !=NULL && inputBufferSize >0 );
    assert( outputBuffer!=NULL && outputBufferSize>0 );
    
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
    
    return Inf_decompress(inflater, (Byte*)outputBuffer, outputBufferSize);
}


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
    
    
    if ( inf->action==InfAction_Finish ) { return 0; }

    while ( destBytesStillNeeded>0 && inf->action!=InfAction_Finish )
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
                inf->providedData=inf->helperInput;
                inf->dataProviderFunc(inf, &inf->providedData );
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
