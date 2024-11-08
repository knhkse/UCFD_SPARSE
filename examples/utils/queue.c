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

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"


/**
 * @brief       Initialize queue structure with NULL value
 * @param       queue       Queue data structure to initialize
 */
void initQueue(struct Queue *queue)
{
    queue -> front = NULL;
    queue -> rear = NULL;
}


/**
 * @brief       Check if queue structure is empty
 * @param       queue       Queue data structure to check
 */
int isEmpty(struct Queue *queue)
{
    if (queue->front == NULL)
        return 1;
    else
        return 0;
}


/**
 * @brief       Allocate new node and add to queue structure
 * @param       queue       Queue data structure to add new node
 * @param       data        Value to add
 */
void enqueue(struct Queue *queue, int data)
{
    struct Node *newNode = (struct Node *)malloc(sizeof(struct Node));
    newNode -> data = data;
    newNode -> next = NULL;
    
    if (isEmpty(queue)) {
        queue -> front = newNode;
        queue -> rear = newNode;
    }
    else {
        queue -> rear -> next = newNode;
        queue -> rear = newNode;
    }
}


/**
 * @brief       Take out the first added element
 * @param       queue       Queue data structure to take element out
 * @param[out]  popdata     Value of poped data
 */
int dequeue(struct Queue *queue)
{
    struct Node *delNode;
    int popdata;
    
    // Raise exception if queue is empty
    if (isEmpty(queue)){
        printf("Empty Queue\n");
        return -1;
    }

    // 
    delNode = queue -> front;
    popdata = delNode -> data;
    queue -> front = queue -> front -> next;

    free(delNode);

    return popdata;
}