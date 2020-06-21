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

namespace PDZip {
    class ReversedHuffmanDecoder;
}

#undef  Byte
#define Byte  unsigned char
#define InfMinBufferSize    (32*1024)
#define InfHelperBufferSize (64*1024)

#define Inf_EndOfBlock           256
#define Inf_MaxValidLengthCode   285
#define Inf_MaxValidDistanceCode 29
#define Inf_LastValidLength      18
#define Inf_LastValidCode        290
#define Inf_CodeLengthTableSize  ((Inf_LastValidLength+1)+(Inf_LastValidCode+1))
#define Inf_NextIndexMask        0x03FF



typedef enum InfError {
    InfError_BadBlockContent = -8,     ///< The content of block is invalid
    InfError_BadBlockLength = -7,      ///< The length of block is invalid
    InfError_BadBlockType = -6,        ///< The type of block is invalid
    InfError_UnsupportedZLibHeader = -5,
    InfError_BadZLibHeader = -4,
    InfError_BadParameter = -3,
    InfError_Adler32Mismatch = -2,
    InfError_Failed = -1,
    InfError_None = 0
} InfStatus;

typedef enum InfAction {
    InfAction_Finish                 = 0,
    InfAction_FillInputBuffer        = 1,
    InfAction_UseOutputBufferContent = 2,
    InfAction_ProcessNextChunk       = 256,
    InfAction_Init                   = 1024
} InfAction;

typedef int InfBool; ///< Boolean value
#define Inf_FALSE 0
#define Inf_TRUE  1


typedef struct InfData {
    const void* buffer;
    size_t      bufferSize;
} InfDatabuffer;


typedef void (*InfDataReceiverFunc)(struct Inflater* inflater, const unsigned char* bytes, size_t numberOfBytes);
typedef void (*InfDataProviderFunc)(struct Inflater* inflater, InfData* data);


typedef struct Inflater {
    /* PDZip::Inflater* obj; */
    
    int       mode;
    int       flags;
    InfAction action;
    InfError  error;
    
    int finished; /* << HACK!!!  */
    
    void*               userPtr;
    InfDataProviderFunc dataProviderFunc;
    InfDataReceiverFunc dataReceiverFunc;
    
    /* inflaterProcessChunk(..) */

    Byte*       outputChunk;             ///< pointer to the last decompressed chunk of data
    size_t      outputChunkSize;         ///< number of decompressed bytes in 'outputChunk'
    size_t      outputBufferContentSize; ///< total number of decompressed bytes in OutputBuffer when 'InfAction_UseOutputBufferContent'

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
    int       step;
    
    unsigned _zlibheader;
    unsigned _zlibheader_method; // < compression method
    unsigned _zlibheader_wsize;  // < window size
    unsigned _zlibheader_level;  // < compression level
    
    unsigned _lastBlock;
    unsigned _blocktype;
    unsigned bitbuffer;  ///< bit-stream buffer
    unsigned bitcounter; ///< number of bits contained in 'bitbuffer'
    
    unsigned _hlit;
    unsigned _hdist;
    unsigned _hclen;
    
    unsigned _literal;
    unsigned _seq_dist;
    unsigned _seq_len;
    
    /* HIDDEN: code lengths reader */
    struct {
        unsigned  command;
        unsigned  code;
        unsigned  length;
        unsigned  repetitions;
        unsigned* insertPtr[Inf_LastValidLength+1];
        unsigned  table[Inf_CodeLengthTableSize];
        unsigned char lengths[19];
        int       nextIndex;
        int       size;
    } cl;


    
    PDZip::ReversedHuffmanDecoder* _frontDecoder;     ///< The base decoder used to decode the next 2 decoders (it's crazy!)
    PDZip::ReversedHuffmanDecoder* _literalDecoder;   ///< The literal+length huffman decoder
    PDZip::ReversedHuffmanDecoder* _distanceDecoder;  ///< The distance huffman decoder


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
