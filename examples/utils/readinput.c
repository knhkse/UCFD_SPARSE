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







