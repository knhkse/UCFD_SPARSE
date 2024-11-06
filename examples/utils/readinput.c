/** ======================================================================================================================
 * @file        queue.c
 * @brief       
 * @details     Queue data structure for computing Reverse Cuthill-McKee algorithm.
 *              Queue is a FIFO (First-In First-Out) data structure to access data in serial input order.  
 * 
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 * 
 * @date        July 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 * 
 * =======================================================================================================================
 */

#include <regex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "readinput.h"

void input_params(FILE *fp, int *vars)
{
    /* File read variables */
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    /* regex variables */
    regex_t pattern;
    regmatch_t groupArray;
    char *scopy = NULL;
    int nwidth;
    int index = 0;

    regcomp(&pattern, "[[:digit:]]+", REG_EXTENDED);

    while((read = getline(&line, &len, fp)) != -1) {
        
        if (regexec(&pattern, line, 1, &groupArray, 0) == 0) {
            scopy = line;
            while (regexec(&pattern, scopy, 1, &groupArray, 0) != REG_NOMATCH) {
                regexec(&pattern, scopy, 1, &groupArray, 0);
                scopy += groupArray.rm_so;
                nwidth = groupArray.rm_eo - groupArray.rm_so;
                char strnum[nwidth + 1];
                strncpy(strnum, scopy, nwidth);
                vars[index] = atoi(strnum);
                scopy += groupArray.rm_eo;
                index++;
            }
        }
    }
    regfree(&pattern);
}







