/*
    main.c
    uzo-proc
  
    Created by Martin on 23/06/2020.
    Copyright © 2020 Martin Rizzo. All rights reserved.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#define INFLATER_IMPLEMENTATION
#include "inflater.h"


#define CHUNK_SIZE 512


#define isSequenceEqual(ptr1, ptr2, numberOfBytes) ( ptr1[0]==ptr2[0] && 0==memcmp(ptr1,ptr2,numberOfBytes) )

static long getFileSize(FILE* file) { fseek(file,0,SEEK_END);  return ftell(file); }

static long min(long a, long b) { return a<b ? a : b; }

/**
 * Finds the position of a sequence of bytes in the provided file
 * @param file                   The file to search
 * @param numberOfBytesToSearch  The maximum number of bytes to search in the provided file,
 * @param sequenceToFind         Pointer to the sequence of bytes to find
 * @param sequenceSize           The number of bytes contained in 'sequenceToFind'
 *                               (-1 = search all the file content)
 * @returns
 *    The position in the file where the sequence was found or (-1) if it was not found
 */
long reverseFindInFile(FILE* file, size_t numberOfBytesToSearch, const char* sequenceToFind, size_t sequenceSize) {
    char buffer[CHUNK_SIZE]; const char* ptr;
    long filePosition, fileSizeLeft; int bytesToRead, bytesToMove;
    const long fileSize        = file!=NULL ? getFileSize(file) : 0;
    const long fileSizeMinimum = numberOfBytesToSearch>=0 && fileSize>numberOfBytesToSearch ? (fileSize-numberOfBytesToSearch) : 0;
    assert( file!=NULL && sequenceToFind!=NULL );
    assert( 1<=sequenceSize && sequenceSize<CHUNK_SIZE );
    
    bytesToMove=0;
    for ( fileSizeLeft=fileSize; fileSizeLeft>fileSizeMinimum; fileSizeLeft=filePosition) {
        /* read a chunk of data from file */
        bytesToRead  = (int)min( (CHUNK_SIZE-bytesToMove) , fileSizeLeft );
        filePosition = (fileSizeLeft - bytesToRead);
        if (bytesToMove>0) { memmove(&buffer[bytesToRead], &buffer[0], bytesToMove); }
        fseek(file, filePosition, SEEK_SET); fread(buffer, 1, bytesToRead, file);
        /* try to find the sequence in the current chunk of data */
        for ( ptr=&buffer[bytesToRead+bytesToMove-sequenceSize]; ptr>=buffer; --ptr) {
            if ( isSequenceEqual(ptr,sequenceToFind,sequenceSize) ) { return filePosition+(ptr-buffer); }
        }
        bytesToMove = (int)sequenceSize;
    }
    return -1;
}


void getFirstCompressedFile(FILE* zipFile) {
    long position;
    static const char eocdSignature[4] = { 0x50, 0x4B, 0x05, 0x06 };
    position = reverseFindInFile(zipFile,65536, eocdSignature,sizeof(eocdSignature));
    printf("EndOfDirectory position = %ld\n", position);
}

void unzipFirstFile(const char* zipFilePath) {
    FILE* zipFile;
    assert( zipFilePath!=NULL && zipFilePath[0]!='\0' );
    
    zipFile = fopen(zipFilePath, "rb");
    if (zipFile==NULL) { return; }
    
    getFirstCompressedFile(zipFile);
        
    fclose(zipFile);
}

/*=================================================================================================================*/
#pragma mark - > MAIN

#define VERSION   "0.1"
#define COPYRIGHT "Copyright (c) 2020 Martin Rizzo"
#define isOption(param,name1,name2) \
    (strcmp(param,name1)==0 || strcmp(param,name2)==0)

/**
 * Application starting point
 * @param argc  The number of elements in the 'argv' array
 * @param argv  An array containing each command-line parameter (starting at argv[1])
 */
int main(int argc, char *argv[]) {
    const char **filePaths; int numberOfFiles;
    const char *param; int i;
    int printHelpAndExit=0, printVersionAndExit=0;
    const char *help[] = {
        "USAGE: uzo-proc [options] file1.zip file2.zip ...","",
        "  OPTIONS:",
        "    -h, --help             display this help and exit",
        "    -v, --version          output version information and exit",
        NULL
    };
    

    /* process all flags/options and store all filenames */
    filePaths     = malloc(argc * sizeof(char*));
    numberOfFiles = 0;
    for (i=1; i<argc; ++i) { param=argv[i];
        if ( param[0]!='-' ) { filePaths[numberOfFiles++]=param; }
        else {
            if      ( isOption(param,"-h","--help")    ) { printHelpAndExit=1;    }
            else if ( isOption(param,"-v","--version") ) { printVersionAndExit=1; }
        }
    }
    
    /* print help or version if requested */
    if ( printHelpAndExit    ) { i=0; while (help[i]!=NULL) { printf("%s\n",help[i++]); } return 0; }
    if ( printVersionAndExit ) { printf("UZO (unzip one file) version %s\n%s\n", VERSION, COPYRIGHT);    return 0; }
    
    /* decompress the first file of each requested zip file */
    for (i=0; i<numberOfFiles; ++i) {
        unzipFirstFile(filePaths[i]);
    }
    
    free(filePaths);
    return 0;
}
