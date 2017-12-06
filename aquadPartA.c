/*
  The first note in my implementation strategy is that initially I break the
  function into numprocs - 1 segments and send each of them to a worker.
  This is opposed to just putting the intitial boundaries on the stack
  and letting the work grow until all the workers have a task. In this way all
  of the workers have work to do from the beginning. At this point the farmer
  enters a loop that ends when the stack is empty and no workers are currently
  working. For communication between the farmer and the worker only MPI_Send
  and MPI_Recv are used.

  The farmer is always receiving messages from ANY_SOURCE and ANY_TAG and then
  uses the source field to determine who to send the next message to. The tag
  field the farmer receives is used to determine whether or not this message is
  the end of a workload or a message to add more work to the bag. If this
  message contains more work the farmer adds 2 new tasks to the bag and
  notes that this worker is active. I noticed that it can happen that a worker
  will return a result (as opposed to more work) and there isn't anything else
  in the stack. In this case the worker will become idle but other workers may
  return more work. So, at the point where the farmer is receiving new work it
  will check to see if there are any processes that are idle and if so give one
  of them one task and the worker that returned the tasks the other task.
  If the worker
  returns a value, and no more work, the farmer adds the value in
  the message to the result and indicates this worker is currently idle. At this
  point the farmer will check to see if there is any more work in the bag and if
  so will send a new job to the worker it last received a message from and note
  that this worker isn't idle. Once all the workers are idle and there are no
  more tasks in the bag the farmer will send a signal to all the workers that
  they are finished.

  The workers very simply wait to receive a message from the farmer. They will
  continue receiving messages from the farmer until they recieve the quit tag.
  Once a worker receives a message it will decompose the message into the
  upper and lower bounds of the function and do its computation. If the result
  is above the specified threshold the worker will send a message back to the
  farmer that contains information to add 2 new tasks to the bag. Otherwise, the
  worker will send a message to the farmer that contains the result obtained.
  Whether this message is new tasks or a result is specified using the tag field
  of the MPI_Send.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "stack.h"

#define EPSILON 1e-3
#define F(arg)  cosh(arg)*cosh(arg)*cosh(arg)*cosh(arg)
#define A 0.0
#define B 5.0

#define FARMER 0
#define ADD_TASK 0
#define FINISHED 1
#define SLEEPTIME 1
#define NO_MORE_WORK 2

int LAST_WORKER;
int *tasks_per_process;

double farmer(int);
void worker(int);
int isdone(double*, int);
void nextidle(int, double*, int);

int main(int argc, char **argv ) {
  int i, myid, numprocs;
  double area, a, b;

  LAST_WORKER = 0;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if(numprocs < 2) {
    fprintf(stderr, "ERROR: Must have at least 2 processes to run\n");
    MPI_Finalize();
    exit(1);
  }

  if (myid == 0) { // Farmer
    // init counters
    tasks_per_process = (int *) malloc(sizeof(int)*(numprocs));
    for (i=0; i<numprocs; i++) {
      tasks_per_process[i]=0;
    }
  }

  if (myid == 0) { // Farmer
    area = farmer(numprocs);
  } else { //Workers
    worker(myid);
  }

  if(myid == 0) {
    fprintf(stdout, "Area=%lf\n", area);
    fprintf(stdout, "\nTasks Per Process\n");
    for (i=0; i<numprocs; i++) {
      fprintf(stdout, "%d\t", i);
    }
    fprintf(stdout, "\n");
    for (i=0; i<numprocs; i++) {
      fprintf(stdout, "%d\t", tasks_per_process[i]);
    }
    fprintf(stdout, "\n");
    free(tasks_per_process);
  }
  MPI_Finalize();
  return 0;
}

double farmer(int numprocs) {
  int tag, i, who;
  double points[2], temp[3], result, step, tempResult, *data, iscomplete[numprocs-1];

  stack *stack = new_stack();
  MPI_Status status;

  result = 0;
  // divide the function into numprocs-1 segments 
  // send them to the workers
  step = (B-A) / (numprocs-1);
  points[0] = A;
  for (i=0; i<numprocs-1; i++) {
    iscomplete[i] = EPSILON + 1;
    points[1] = points[0] + step;
    MPI_Send(points, 2, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
    points[0] = points[1];
  }

  // continue receiving messages until we having computed the area
  while (!isdone(iscomplete, numprocs) || !is_empty(stack)) {
    MPI_Recv(temp, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
    who = status.MPI_SOURCE;

    // received new tasks - add them to the bag
    if (tag == ADD_TASK) {
      points[0] = temp[0];
      points[1] = temp[1];
      push(points, stack);
      points[0] = points[1];
      points[1] = temp[2];
      push(points, stack);

      // check for another idle worker
      nextidle(who, iscomplete, numprocs);
      if (LAST_WORKER != who) {
	data = pop(stack);
	iscomplete[LAST_WORKER-1] = EPSILON + 1;
	tasks_per_process[LAST_WORKER] += 1;
	MPI_Send(data, 2, MPI_DOUBLE, LAST_WORKER, 0, MPI_COMM_WORLD);
      }
      iscomplete[who-1] = EPSILON + 1;
    }
    // no new task - add to the result
    else {
      result += temp[0];
      iscomplete[who-1] = EPSILON - 1;
    }
    // if there are any tasks send one to this worker
    if (!is_empty(stack)) {
      data = pop(stack);
      iscomplete[who-1] = EPSILON + 1;
      tasks_per_process[status.MPI_SOURCE] += 1;
      MPI_Send(data, 2, MPI_DOUBLE, who, 0, MPI_COMM_WORLD);
    }
  }
  // tell all the workers to finish
  for (i=0; i<numprocs-1; i++)
    MPI_Send(data, 2, MPI_DOUBLE, i+1, NO_MORE_WORK, MPI_COMM_WORLD);

  return result;
}

void worker(int mypid) {
  int tag;
  MPI_Status status;
  double mid, fmid, larea, rarea, lrarea, left, right, fleft, fright, task[3];

  MPI_Recv(task, 2, MPI_DOUBLE, FARMER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  tag = status.MPI_TAG;
  while (tag != NO_MORE_WORK) {
    usleep(SLEEPTIME);

    // get left, right, and mid bounds
    // compute the function value at those points
    left = task[0]; right = task[1];
    fleft = F(left); fright = F(right);
    mid = (left + right) / 2;
    fmid = F(mid);

    // this is the previous area
    lrarea = (fleft+fright) * (right-left)/2;

    // get the areas of the left half and right half
    larea = (fleft + fmid) * (mid - left) / 2;
    rarea = (fmid + fright) * (right - mid) / 2;

    // result is above threshold - send back new tasks
    if( fabs((larea + rarea) - lrarea) > EPSILON ) {
      task[0] = left;
      task[1] = mid;
      task[2] = right;

      MPI_Send(task, 3, MPI_DOUBLE, FARMER, ADD_TASK, MPI_COMM_WORLD);
    }
    // send back result
    else {
      task[0] = larea + rarea;
      MPI_Send(task, 1, MPI_DOUBLE, FARMER, FINISHED, MPI_COMM_WORLD);
    }
    // wait for next message
    MPI_Recv(task, 2, MPI_DOUBLE, FARMER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
  }
}

// takes the array of idle workers and number of processes
// returns whether or not all the processes are idle
int isdone(double *iscomplete, int numprocs) {
  int done = 1, i;
  for (i=0; i<numprocs-1; i++)
      if (iscomplete[i] > EPSILON)
	done = 0;

  return done;
}

// takes the start location in the array of idle processes, the array of idle
// processes, and the total number of processes
// sets LAST_WORKER to the first worker it finds idle after the start 
void nextidle(int start, double *iscomplete, int numprocs) {
  int i;
  for (i=start-1; i<(numprocs-1)+start-1; i++) {
    if (iscomplete[i%(numprocs-1)] < EPSILON) {
      LAST_WORKER = i;
      break;
    }
  }
  LAST_WORKER = i%(numprocs-1)+1;
}
