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


#define Inf_LastValidLength 18
#define Inf_InvalidLength   24
#define Inf_MaxNumberOfCodes 19
#define Inf_CodeLengthTableSize ((Inf_LastValidLength+1)+(Inf_LastValidCode+1))

#define inf (*infptr)

/* Macros used for control flow inside the 'inflateProcessChunk(..) function */
#define inf__FILL_INPUT_BUFFER()                             inf.action=InfAction_FillInputBuffer;        break;
#define inf__USE_OUTPUT_BUFFER_CONTENT()                     inf.action=InfAction_UseOutputBufferContent; break;
#define inf__FINISH()               step=InfStep_End;        inf.action=InfAction_Finish;                 break;
#define inf__ERROR(err)             step=InfStep_FatalError; inf.action=InfAction_Finish; inf.error=err;  break;
#define inf__goto(dest_step)        step=dest_step;                                                       break;
#define inf__fallthrough(next_step) step=next_step;


/*====================================================================================================================*/
#   pragma mark - REVERSE BITS

static const unsigned char Inf_Reverse[256] = {
  0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};



/*====================================================================================================================*/
#   pragma mark - HUFFMAN DECODER


#define InfHuff_ReverseBits(huffman, length) (                                                    \
  ((InfHuff_ReverseArray[huffman&0xFF]<<8) | (InfHuff_ReverseArray[huffman>>8])) >> (16-(length)) \
)

#define InfHuff_NextCanonical(huffman, currentLength, newLength)  \
    huffman       = (huffman+1) << (newLength - currentLength);   \
    huffmanLength = newLength                                     \


/**
 * Fills a huffman table (or subtable) with all posible values that match the provided huffman code
 * @param table            Pointer to the huffman table (or sub-table) to fill
 * @param tableSize        Number of entries available in `table`
 * @param huffman          The huffman code
 * @param huffmanLength    The length of the huffman code in number of bits
 * @param byte             The decoded byte corresponding to the provided huffman code
 * @param discardedBits    Number of bits to discard (0=table, 8=sub-table)
 */
static void InfHuff_FillTable(InfHuff*  table,
                              unsigned  tableSize,
                              unsigned  huffman,
                              unsigned  huffmanLength,
                              unsigned  byte,
                              unsigned  discardedBits)
{
    unsigned unknownBits, knownHuffman; InfHuff data;
    const unsigned unknownStep    = (1<<(huffmanLength-discardedBits));
    const unsigned reverseHuffman = ((Inf_Reverse[huffman&0xFF]<<8)|(Inf_Reverse[huffman>>8])) >> (16-huffmanLength);
    assert( table!=NULL && huffmanLength>discardedBits && reverseHuffman<(1<<huffmanLength) );
    
    InfHuff_Set(data, reverseHuffman, huffmanLength, byte);
    knownHuffman = reverseHuffman >> discardedBits;
    for (unknownBits=0; unknownBits<tableSize; unknownBits+=unknownStep) {
        table[ unknownBits | knownHuffman ] = data;
    }
}

/**
 * Makes a huffman table to be used in fast extraction from bitstream
 * @param table    Pointer to the huffman table to fill
 * @param codelen  The first element of a list containing `code-length` data sorted by length
 * @returns
 *      The same table pointer provided in the first parameter.
 */
const InfHuff* InfHuff_MakeTable(InfHuff* table, const InfCodelen *codelen) {
    static const InfHuff invalid = InfHuff_Const(0, Inf_InvalidLength);
    InfHuff data; int i, subtableIndex; const InfCodelen *codelen_end;
    unsigned subtableSize, huffman, huffmanLength;
    assert( table!=NULL && codelen!=NULL  );
    
    /* reset the main-table */
    for (i=0; i<Inf_MainTableSize; ++i) { table[i] = invalid; }
    
    /* init huffman canonical code */
    huffman       = 0;
    huffmanLength = codelen->length;

    /* lengths from 1 to 8                              */
    /* unknown bits are filled with all possible values */
    while ( codelen!=NULL && huffmanLength<=8 ) {
        InfHuff_FillTable(table, 256, huffman, huffmanLength, codelen->code, 0);
        if ( (codelen=codelen->next)!=NULL ) { InfHuff_NextCanonical(huffman, huffmanLength, codelen->length); }
    }
    /* lengths from 9         */
    /* subtables are created  */
    subtableIndex = Inf_MainTableSize;
    while ( codelen!=NULL )
    {
        const InfCodelen* const codelen_first = codelen;
        const unsigned          huffman_first = huffman;
        const unsigned          index         = huffman >> (huffmanLength-8);
        
        /* calculate subtable size (and find first and last element) */
        do {
            subtableSize = huffmanLength;
            if ( (codelen=codelen->next)!=NULL ) { InfHuff_NextCanonical(huffman, huffmanLength, codelen->length); }
        } while ( codelen!=NULL && (huffman>>(huffmanLength-8))==index );
        codelen_end  = codelen;
        subtableSize = (1<<(subtableSize-8));
        assert( subtableIndex+subtableSize <= Inf_HuffTableSize );
        
        /* create subtable */
        huffman       = huffman_first;
        huffmanLength = codelen_first->length;
        table[Inf_Reverse[index]] = InfHuff_SubTableRef(data, subtableIndex, subtableSize-1);
        codelen = codelen_first; while ( codelen!=codelen_end ) {
            InfHuff_FillTable(&table[subtableIndex], subtableSize, huffman, huffmanLength, codelen->code, 8);
            if ( (codelen=codelen->next)!=NULL ) { InfHuff_NextCanonical(huffman, huffmanLength, codelen->length); }
        }
        subtableIndex += subtableSize;
    }
    return table;
}


#define InfHuff_Decode(table, tmp, bitbuffer) \
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
    data=InfHuff_Decode(table, tmp, inf.bitbuffer); numberOfBitsToAdvance=data.value.length;
    
    while (inf.bitcounter<numberOfBitsToAdvance) {
        if (!InfBS_LoadNextByte(infptr)) { return Inf_FALSE; }
        data=InfHuff_Decode(table, tmp, inf.bitbuffer); numberOfBitsToAdvance=data.value.length;
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
    (*dest) = inf.bitbuffer;
    inf.bitbuffer  = inf.bitcounter = 0;
    return Inf_TRUE;
}

/** Reads a sequence of bytes directly from the input stream (bitbuffer must be empty) */
static InfBool InfBS_ReadBytes(Inflater* infptr, unsigned char* dest, size_t* inout_numberOfBytes) {
    size_t numberOfRequestedBytes, successfulBytes;
    assert( infptr!=NULL && dest!=NULL && inout_numberOfBytes!=NULL && inf.bitcounter==0 );
    numberOfRequestedBytes = (*inout_numberOfBytes);
    successfulBytes        = min( numberOfRequestedBytes, (inf.inputChunkEnd-inf.inputChunkPtr) );
    memcpy(dest, inf.inputChunkPtr, successfulBytes);
    inf.inputChunkPtr += successfulBytes;
    (*inout_numberOfBytes) = successfulBytes;
    return (successfulBytes == numberOfRequestedBytes);
}


/*====================================================================================================================*/
#   pragma mark - READING CODE/LENGTHS PAIRS


static void InfCL_Add(Inflater* infptr, int code, unsigned length) {
    assert( 0<=code   &&   code<=Inf_LastValidCode   );
    assert( 0<=length && length<=Inf_LastValidLength );
    if (length>0) {
        InfCodelen* const newElement = &inf.clList.elements[ inf.clList.index++ ];
        InfCodelen* const last       = inf.clList.tailPtr[length];
        if (last) { last->next = newElement; } else { inf.clList.headPtr[length] = newElement; }
        newElement->code   = code;
        newElement->length = length;
        newElement->next   = NULL;
        inf.clList.tailPtr[length] = newElement;
    }
}

static void InfCL_AddRange(Inflater* infptr, int firstCode, int lastCode, unsigned length) {
    int code;
    assert( infptr!=NULL );
    assert( 0<=firstCode && firstCode<=(Inf_LastValidCode)   );
    assert( 0<=lastCode  &&  lastCode<=(Inf_LastValidCode+1) );
    assert( firstCode<=lastCode );
    for (code=firstCode; code<lastCode; ++code) {
        InfCL_Add(infptr, code,length);
    }
}


static void InfCL_Open(Inflater* infptr, InfBool resetRepetitions) {
    int length;
    if (resetRepetitions) { inf.cl.command = inf.cl.length = inf.cl.repetitions = 0; }
    inf.cl.code      = 0;
    inf.clList.index = 0;
    for (length=0; length<=Inf_LastValidLength; ++length) {
        inf.clList.headPtr[length] = inf.clList.tailPtr[length] = NULL;
    }
}

static const InfCodelen* InfCL_Close(Inflater* infptr) {
    int lastValidLength = 0;
    int currentLength   = 0;
    InfCodelen *firstElement, *currentHead;
    
    /* connect all lists creating a big one sorted by `length` */
    while (inf.clList.headPtr[currentLength]==NULL) { ++currentLength; }
    firstElement = inf.clList.headPtr[ lastValidLength = currentLength ];
    
    do {
        do { currentHead = inf.clList.headPtr[ ++currentLength ];
        } while ( currentLength<Inf_LastValidLength && currentHead==NULL );
        inf.clList.tailPtr[lastValidLength]->next = currentHead;
        lastValidLength = currentLength;
    } while ( currentLength<Inf_LastValidLength );
    inf.clList.index = 0;
    return firstElement;
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
    for (i=0; i<Inf_MaxNumberOfCodes; i++) { InfCL_Add(infptr, i, inf.cl.lengths[i]);  }
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
        InfCL_Add(infptr, inf.cl.code, inf.cl.length);
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
            InfCL_Add(infptr, inf.cl.code, inf.cl.length);
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
                InfCL_Add(infptr, inf.cl.code, inf.cl.length);
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


static InfHuff* getFixedLiteralDecoder(Inflater* infptr) {
    static InfHuff table[Inf_HuffTableSize];
    static InfBool loaded = Inf_FALSE;
    if (!loaded) {
        loaded = Inf_TRUE;
        InfCL_Open(infptr, Inf_TRUE);
        InfCL_AddRange(infptr,   0,144, 8);
        InfCL_AddRange(infptr, 144,256, 9);
        InfCL_AddRange(infptr, 256,280, 7);
        InfCL_AddRange(infptr, 280,288, 8);
        InfHuff_MakeTable(table, InfCL_Close(infptr) );
    }
    return table;
}

static InfHuff* getFixedDistanceDecoder(Inflater* infptr) {
    static InfHuff table[Inf_HuffTableSize];
    static InfBool loaded = Inf_FALSE;
    if (!loaded) {
        loaded = Inf_TRUE;
        InfCL_Open(infptr, Inf_TRUE);
        InfCL_AddRange(infptr, 0,32, 5);
        InfHuff_MakeTable(table, InfCL_Close(infptr) );
    }
    return table;
}





/*====================================================================================================================*/
#   pragma mark - IMPLEMENTATION OF PUBLIC FUNCTIONS


Inflater* inflaterCreate(void* workingBuffer, size_t workingBufferSize) {
    Inflater* infptr = (Inflater*)malloc(sizeof(Inflater));
    
    inf.mode            = InfMode_Uninitialized;
    inf.flags           = InfFlags_FreeItself;
    inf.action          = InfAction_Init;
    inf.error           = InfError_None;
    
    inf.userPtr          = NULL;
    inf.dataProviderFunc = NULL;
    inf.dataReceiverFunc = NULL;
    
    inf.inputChunkPtr   = NULL;
    inf.inputChunkEnd   = NULL;
    inf.outputChunk     = NULL;
    inf.outputChunkSize = 0;
    inf.outputBufferContentSize = 0;
    
    /*---- inflaterTake / inflaterFeed ----------- */
    inf.helperOutput.buffer     = NULL;
    inf.helperOutput.bufferSize = 0;
    inf.helperInput.buffer      = NULL;
    inf.helperInput.bufferSize  = 0;
    inf.takeOutputPtr           = NULL;
    inf.takeOutputRemaining     = 0;
    inf.providedData.buffer     = NULL;
    inf.providedData.bufferSize = 0;
    
    /*---- decompress ---------------------------- */
    inf.step       = 0;
    inf.bitbuffer  = 0;
    inf.bitcounter = 0;
    
    
    if (workingBuffer!=NULL && workingBufferSize>=InfHelperBufferSize) {
        inf.helperOutput.buffer     = workingBuffer;
        inf.helperOutput.bufferSize = workingBufferSize;
    }
    return infptr;
}

/**
 * Decompresses a chunk of compressed data
 * @param infptr            Pointer to the `Inflater` object created with `inflaterCreate(..)`
 * @param outputBuffer      Pointer to the destination buffer where decompressed data will be stored
 * @param outputBufferSize  The available capacity of `destBuffer` in number of bytes
 * @param inputBuffer       Pointer to the source buffer from where compressed data is read
 * @param inputBufferSize   The length of `sourBuffer` in number of bytes
 * @returns
 *     The action required to execute to continue with the next chunk, ex: `InfAction_FillInputBuffer`
 */
InfAction inflaterProcessChunk(Inflater*         infptr,
                               void* const       outputBuffer,
                               size_t            outputBufferSize,
                               const void* const inputBuffer,
                               size_t            inputBufferSize)
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
    
    InfBool canReadAll; InfStep step; size_t numberOfBytes; unsigned temp; unsigned char *sequencePtr, *writePtr;
    unsigned char* const writeEnd = ((Byte*)outputBuffer + outputBufferSize);
    
    assert( inputBuffer !=NULL && inputBufferSize >0 );
    assert( outputBuffer!=NULL && outputBufferSize>0 );

    switch ( inf.action ) {
        case InfAction_Init:
            inf.outputChunk   = (Byte      *)outputBuffer;
            inf.inputChunkPtr = (const Byte*)inputBuffer;
            inf.inputChunkEnd = inf.inputChunkPtr + inputBufferSize;
            break;
        case InfAction_FillInputBuffer:
            inf.outputChunk  += inf.outputChunkSize;
            inf.inputChunkPtr = (const Byte*)inputBuffer;
            inf.inputChunkEnd = inf.inputChunkPtr + inputBufferSize;
            break;
        case InfAction_UseOutputBufferContent:
            inf.outputChunk = (Byte*)outputBuffer;
            inf.outputBufferContentSize = 0;
            break;
        case InfAction_ProcessNextChunk:
            inf.outputChunk += inf.outputChunkSize;
            break;
        default: /* InfAction_Finish */
            inf.outputChunkSize = 0;
            return inf.action;
    }

    assert( inf.inputChunkPtr!=NULL );
    assert( inf.inputChunkPtr<=inf.inputChunkEnd );
    assert( inf.outputChunk!=NULL && inf.outputChunk<writeEnd );
    assert( outputBuffer!=NULL );
    assert( outputBufferSize>=(32*1024) );
    
    step       = (InfStep)inf.step;
    writePtr   = inf.outputChunk;
    inf.action = InfAction_ProcessNextChunk;
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
                inf._seq_len = (temp & 0xFFFF);
                if ( inf._seq_len != ((~temp)>>16) ) { inf__ERROR(InfError_BadBlockLength); }
                inf__fallthrough(InfStep_Output_UncompressedBlock);
                
            case InfStep_Output_UncompressedBlock:
                numberOfBytes = min( inf._seq_len, (writeEnd-writePtr) );
                canReadAll    = InfBS_ReadBytes(infptr, writePtr, &numberOfBytes);
                inf._seq_len -= numberOfBytes;
                writePtr     += numberOfBytes;
                if      ( !canReadAll    ) { inf__FILL_INPUT_BUFFER();         }
                else if ( inf._seq_len>0 ) { inf__USE_OUTPUT_BUFFER_CONTENT(); }
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
                inf.frontDecoder = InfHuff_MakeTable(inf.huffmanTable1, InfCL_Close(infptr));
                inf__fallthrough(InfStep_Load_LiteralHuffmanTable);
                
            /*----------------------------------
             * load literal-length huffman table
             */
            case InfStep_Load_LiteralHuffmanTable:
                InfCL_Open(infptr,Inf_TRUE);
                inf__fallthrough(InfStep_Load_LiteralHuffmanTable2);
                
            case InfStep_Load_LiteralHuffmanTable2:
                if ( !InfCL_ReadCompressedCodes(infptr,inf._hlit+257,inf.frontDecoder) ) { inf__FILL_INPUT_BUFFER(); }
                inf.literalDecoder = InfHuff_MakeTable(inf.huffmanTable0, InfCL_Close(infptr));
                inf__fallthrough(InfStep_Load_DistanceHuffmanTable);
                
            /*----------------------------------
             * load distance huffman table
             */
            case InfStep_Load_DistanceHuffmanTable:
                InfCL_Open(infptr,Inf_FALSE);
                inf__fallthrough(InfStep_Load_DistanceHuffmanTable2);
                
            case InfStep_Load_DistanceHuffmanTable2:
                if ( !InfCL_ReadCompressedCodes(infptr,inf._hdist+1,inf.frontDecoder) ) { inf__FILL_INPUT_BUFFER(); }
                inf.distanceDecoder = InfHuff_MakeTable(inf.huffmanTable1, InfCL_Close(infptr));
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
                if (sequencePtr<(Byte*)outputBuffer) {
                    sequencePtr += outputBufferSize;
                    /*-- copy bytes ---*/
                    inf._seq_len -= numberOfBytes = min( inf._seq_len, (writeEnd-sequencePtr) );
                    memmove( writePtr, sequencePtr, numberOfBytes ); writePtr+=numberOfBytes;
                    /*-----------------*/
                    if ( inf._seq_len==0    ) { inf__goto(InfStep_Read_LiteralOrLength); }
                    if ( writePtr==writeEnd ) { inf__USE_OUTPUT_BUFFER_CONTENT();        }
                    sequencePtr = outputBuffer;
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


/**
 * Decompresses a chunk of data
 * @param infptr                Pointer to the `Inflater` object created with `inflaterCreate(..)`
 * @param decompressedData      Pointer to the destination buffer where decompressed data will be stored
 * @param decompressedDataSize  The available capacity of `destBuffer` in number of bytes
 */
size_t inflaterTake(Inflater* infptr, void* decompressedData, size_t decompressedDataSize) {
    /* InfAction action; */
    size_t sourAvailableBytes,destBytesStillNeeded;
    Byte* dest; const Byte* sour;
    Byte* outputBuffer; size_t outputBufferSize;
    int fillInputBuffer;

    /* initialization */
    assert( (inf.mode==InfMode_Uninitialized) || (inf.mode==InfMode_Take) );
    if (inf.mode==InfMode_Uninitialized) {
        inf.mode= InfMode_Take;
        if (inf.helperOutput.buffer==NULL) {
            inf.helperOutput.buffer = malloc( inf.helperOutput.bufferSize=InfHelperBufferSize );
            inf.flags |= InfFlags_FreeOutputBuffer;
        }
        if (inf.helperInput.buffer==NULL) {
            inf.helperInput.buffer  = malloc( inf.helperInput.bufferSize=InfHelperBufferSize );
            inf.flags |= InfFlags_FreeInputBuffer;
        }
        inf.takeOutputPtr           = (Byte*)inf.helperOutput.buffer;
        inf.takeOutputRemaining     = 0;
        inf.providedData.buffer     = NULL;
        inf.providedData.bufferSize = 0;
    }

    outputBuffer         = (Byte*)inf.helperOutput.buffer;
    outputBufferSize     =        inf.helperOutput.bufferSize;
    sour                 =        inf.takeOutputPtr;
    sourAvailableBytes   =        inf.takeOutputRemaining;
    dest                 = (Byte*)decompressedData;
    destBytesStillNeeded =        decompressedDataSize;
    
    
    if ( inf.action==InfAction_Finish ) { return 0; }

    while ( destBytesStillNeeded>0 && inf.action!=InfAction_Finish )
    {
        fillInputBuffer = (inf.action==InfAction_Init || inf.action==InfAction_FillInputBuffer);
        
        if (inf.action==InfAction_UseOutputBufferContent && sourAvailableBytes>0 ) {
            /* take available bytes */
            const size_t size = sourAvailableBytes<destBytesStillNeeded ? sourAvailableBytes : destBytesStillNeeded;
            memcpy(dest,sour,size); dest+=size; destBytesStillNeeded-=size; sour+=size; sourAvailableBytes-=size;
            if (sourAvailableBytes==0) { fillInputBuffer=1; sour=outputBuffer; }
        }
        if ( fillInputBuffer ) {
            /* get compressed data from the data provider (when required) */
            if (inf.action==InfAction_Init || inf.action==InfAction_FillInputBuffer) {
                inf.providedData=inf.helperInput;
                inf.dataProviderFunc(infptr, &inf.providedData );
                if ( inf.providedData.bufferSize==0 ) { break; }
            }
            /* decompress the provided data and generate more sourAvailableBytes */
            do {
                inflaterProcessChunk(infptr, outputBuffer, outputBufferSize, inf.providedData.buffer, inf.providedData.bufferSize);
                sourAvailableBytes += inf.outputChunkSize;
            } while ( inf.action==InfAction_ProcessNextChunk );
        }
    }
    inf.takeOutputPtr=sour; inf.takeOutputRemaining=sourAvailableBytes;
    return (decompressedDataSize-destBytesStillNeeded);
}

/**
 * Decompresses a chunk of data
 * @param infptr              Pointer to the `Inflater` object created with `infalterCreate(..)`
 * @param compressedData      Pointer to the source buffer from where compressed data is read
 * @param compressedDataSize  The length of `compressedData` in number of bytes
 */
size_t inflaterFeed(Inflater* infptr, const void* compressedData, size_t compressedDataSize) {
    InfAction action; /* size_t numberOfConsumedBytes=0; */
    Byte* outputBuffer; size_t outputBufferSize;
    assert( (inf.mode==InfMode_Uninitialized) || (inf.mode=InfMode_Feed) );
    assert( inf.dataReceiverFunc!=NULL );

    /* initialization */
    if (inf.mode==InfMode_Uninitialized) {
        inf.mode= InfMode_Feed;
        if (inf.helperOutput.buffer==NULL) {
            inf.helperOutput.buffer = (Byte*)malloc( inf.helperOutput.bufferSize=InfHelperBufferSize );
            inf.flags = InfFlags_FreeOutputBuffer;
        }
    }
    
    /* decompression */
    outputBuffer     = (Byte*)inf.helperOutput.buffer;
    outputBufferSize =        inf.helperOutput.bufferSize;
    do {
        action = inflaterProcessChunk(infptr, outputBuffer, outputBufferSize, compressedData, compressedDataSize);
        if (action==InfAction_UseOutputBufferContent) {
            inf.dataReceiverFunc(infptr, outputBuffer, inf.outputBufferContentSize);
        }
        /* numberOfConsumedBytes += inf.inputChunkSize; */
    } while (action==InfAction_ProcessNextChunk || action==InfAction_UseOutputBufferContent);
    
/*  return numberOfConsumedBytes;  */
    return compressedDataSize;
}

/* receiverFunc     The function that will be called to store the resulting decompressed data */


void inflaterDestroy(Inflater* infptr) {
    if (infptr) {
        /* delete inf.obj; */
        if (0!=(inf.flags & InfFlags_FreeInputBuffer  )) { free((void*)inf.helperInput.buffer);  }
        if (0!=(inf.flags & InfFlags_FreeOutputBuffer )) { free((void*)inf.helperOutput.buffer); }
        if (0!=(inf.flags & InfFlags_FreeItself       )) { free((void*)infptr);                  }
    }
}




#undef inf__FILL_INPUT_BUFFER
#undef inf__USE_OUTPUT_BUFFER_CONTENT
#undef inf__FINISH
#undef inf__ERROR
#undef inf__goto
#undef inf__fallthrough

#undef inf
