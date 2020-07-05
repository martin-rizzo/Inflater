/**
 * @file       inflater.h
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
#ifndef INFLATER_H_INCLUDED
#define INFLATER_H_INCLUDED
#include <stdlib.h>

#undef  Byte
#define Byte  unsigned char
#define InfMinBufferSize    (32*1024)
#define InfHelperBufferSize (64*1024)

#define Inf_EndOfBlock           256
#define Inf_MaxValidLengthCode   285
#define Inf_MaxValidDistanceCode 29
#define Inf_LastValidLength      18
#define Inf_LastValidSymbol      290
#define Inf_CodeLengthMapSize    ((Inf_LastValidLength+1)+(Inf_LastValidSymbol+1))
#define Inf_NextIndexMask        0x03FF


#define Inf_MainTableSize 256       /**< 8bits                        */
#define Inf_HuffTableSize  (2*1024) /**< main-table + all sub-tables  */


typedef enum InfError {
    InfError_BadBlockContent = -8,     /**< The content of block is invalid */
    InfError_BadBlockLength = -7,      /**< The length of block is invalid  */
    InfError_BadBlockType = -6,        /**< The type of block is invalid    */
    InfError_UnsupportedZLibHeader = -5,
    InfError_BadZLibHeader = -4,
    InfError_BadParameter = -3,
    InfError_Adler32Mismatch = -2,
    InfError_Failed = -1,
    InfError_None = 0
} InfError;

typedef enum InfAction {
    InfAction_Finish                 = 0,
    InfAction_FillInputBuffer        = 1,
    InfAction_UseOutputBufferContent = 2,
    InfAction_ProcessNextChunk       = 256,
    InfAction_Init                   = 1024
} InfAction;

typedef int InfBool; /**< Boolean value */
#define Inf_FALSE 0
#define Inf_TRUE  1


typedef struct InfData {
    const void* buffer;
    size_t      bufferSize;
} InfData;


typedef struct InfSymlen {
    unsigned           symbol;        /**< The decompressed value assigned to the huffman code */
    unsigned           huffmanLength; /**< The huffman code length (in number of bits)         */
    struct InfSymlen* next;
} InfSymlen;


typedef union InfHuff {
    struct value    { unsigned length:15, isvalid:1,  code:15; } value;
    struct subtable { unsigned   mask:15,   error:1, index:15; } subtable;
    unsigned raw;
} InfHuff;

#define InfHuff_Set(s, huff, hufflen, byte) s.value.code=byte; s.value.isvalid=1; s.value.length=hufflen
#define InfHuff_SubTableRef(s, index_,mask_) (s.subtable.index=(index_), s.subtable.error=0, s.subtable.mask=(mask_), s)
#define InfHuff_Const(code_,length_) { length_, 1, code_ }

struct Inflater;
typedef void (*InfDataReceiverFunc)(struct Inflater* inflater, const unsigned char* bytes, size_t numberOfBytes);
typedef void (*InfDataProviderFunc)(struct Inflater* inflater, InfData* data);


typedef struct Inflater {

    int       mode;
    int       flags;
    InfAction action;
    InfError  error;
    
    void*               userPtr;
    InfDataProviderFunc dataProviderFunc;
    InfDataReceiverFunc dataReceiverFunc;
    
    /* inflaterProcessChunk(..) */

    Byte*       outputChunk;             /**< pointer to the last decompressed chunk of data  */
    size_t      outputChunkSize;         /**< number of decompressed bytes in 'outputChunk'   */
    size_t      outputBufferContentSize; /**< total number of decompressed bytes in OutputBuffer when 'InfAction_UseOutputBufferContent' */

    const Byte* inputChunkPtr;
    const Byte* inputChunkEnd;
    /* size_t      inputChunkSize; */
    
    
    /* inflaterTake / inflaterFeed */
    
    InfData     helperInput;
    InfData     helperOutput;
    const Byte* takeOutputPtr;
    size_t      takeOutputRemaining;
    InfData     providedData;
    
    
    /* HIDDEN: decompress */
    unsigned step;
    
    unsigned _lastBlock;
    unsigned _blocktype;
    unsigned bitbuffer;  /**< bit-stream buffer                       */
    unsigned bitcounter; /**< number of bits contained in 'bitbuffer' */
    
    unsigned _hlit;
    unsigned _hdist;
    unsigned _hclen;
    
    unsigned _literal;
    unsigned _seq_dist;
    unsigned _seq_len;
    
    /* HIDDEN: symbol-length reader */
    struct {
        InfSymlen*    tailPtr[Inf_LastValidLength+1];
        InfSymlen*    headPtr[Inf_LastValidLength+1];
        InfSymlen     elements[Inf_LastValidSymbol+1]; /**< Elements to add to the list */
        int           elementIndex;        /**< Index to the next free element that is ready to add to the list */
    } symlenList;
    struct {
        unsigned      command;             /**< current command, ex: Command_CopyPreviousLength       */
        unsigned      symbol;              /**< current symbol value                                  */
        unsigned      huffmanLength;       /**< last huffman-length read                              */
        unsigned      repetitions;         /**< number of repetitions of the last huffman-length read */
        unsigned char lengthsBySymbol[19]; /**< Array used to sort lengths by symbol number           */
    } reader;


    const union InfHuff* frontDecoder;        /**< The base decoder used to decode the next 2 decoders (it's crazy!) */
    const union InfHuff* literalDecoder;      /**< The literal+length huffman decoder  */
    const union InfHuff* distanceDecoder;     /**< The distance huffman decoder        */
    InfHuff huffmanTable0[Inf_HuffTableSize];
    InfHuff huffmanTable1[Inf_HuffTableSize];

} Inflater;






extern Inflater* inflaterCreate(void* workingBuffer, size_t workingBufferSize);
extern InfAction inflaterProcessChunk(Inflater* inflater, void* outputBuffer, size_t outputBufferSize, const void* inputBuffer, size_t inputBufferSize);
extern size_t    inflaterTake(Inflater* inflater, void* dest, size_t destSize);
extern size_t    inflaterFeed(Inflater* inflater, const void* sour, size_t sourSize);
extern void      inflaterDestroy(Inflater* inflater);


/*=================================================================================================================*/
#pragma mark - > INTERNAL PRIVATE FUNCTIONS

#ifdef INFLATER_IMPLEMENTATION


#endif /* ifdef INFLATER_IMPLEMENTATION */


#endif /* ifndef INFLATER_H_INCLUDED */
