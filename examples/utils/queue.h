#ifndef QUEUE_H
#define QUEUE_H

/**
 * @brief       Node data structure contains data and next node pointer
 */
struct Node {
    int data;               /** Value */
    struct Node *next;      /** Linked node pointer */
};


/**
 * @brief       Queue data structure
 */
struct Queue {
    struct Node *front;     /** First node element */
    struct Node *rear;      /** Last node element */
};

void initQueue(struct Queue *queue);

int isEmpty(struct Queue *queue);

void enqueue(struct Queue *queue, int data);

int dequeue(struct Queue *queue);


#endif