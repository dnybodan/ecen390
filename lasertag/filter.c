/*
This software is provided for student assignment use in the Department of
Electrical and Computer Engineering, Brigham Young University, Utah, USA.
Users agree to not re-host, or redistribute the software, in source or binary
form, to other persons or other institutions. Users may modify and use the
source code for personal or educational use.
For questions, contact Brad Hutchings or Jeff Goeders, https://ece.byu.edu/
*/

#include "filter.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

// Number of FIR Coeefficitnets
#define FIR_COEF_COUNT 81
//number of queues in z
#define NUM_Z_QUEUES 10
//the z queue size
#define Z_QUEUE_SIZE 10
//default array name size
#define DEFAULT_NAME_SIZE 50
//x queue size
#define X_QUEUE_SIZE 81
//y queue size
#define Y_QUEUE_SIZE 11
//output queue size
#define OUTPUT_QUEUE_SIZE 2000
//decimation factor
#define DECIMATION_FACTOR 10
//iir coefficient count
#define IIR_COEF_COUNT 11
//the power value size
#define POWER_VAL_SIZE 10
//iir a coefficient count
#define IIR_A_COEF_COUNT 10
// Filtering routines for the laser-tag project.
// Filtering is performed by a two-stage filter, as described below.

// 1. First filter is a decimating FIR filter with a configurable number of
// taps and decimation factor.
// 2. The output from the decimating FIR filter is passed through a bank of
// 10 IIR filters. The characteristics of the IIR filter are fixed.

/*********************************************************************************************************
****************************************** Main Filter Functions
******************************************
**********************************************************************************************************/
// declare the queues for filter design
static queue_t xQueue;
static queue_t yQueue;
static queue_t zQueue[FILTER_FREQUENCY_COUNT];
static queue_t outputQueue[FILTER_FREQUENCY_COUNT];
// initialize constants for filter coefficients
const static double fir_coeffs[FIR_COEF_COUNT] = {
    6.0546138291252597e-04,  5.2507143315267811e-04,  3.8449091272701525e-04,
    1.7398667197948182e-04,  -1.1360489934931548e-04, -4.7488111478632532e-04,
    -8.8813878356223768e-04, -1.3082618178394971e-03, -1.6663618496969908e-03,
    -1.8755700366336781e-03, -1.8432363328817916e-03, -1.4884258721727399e-03,
    -7.6225514924622853e-04, 3.3245249132384837e-04,  1.7262548802593762e-03,
    3.2768418720744217e-03,  4.7744814146589041e-03,  5.9606317814670249e-03,
    6.5591485566565593e-03,  6.3172870282586493e-03,  5.0516421324586546e-03,
    2.6926388909554420e-03,  -6.7950808883015244e-04, -4.8141100026888716e-03,
    -9.2899200683230643e-03, -1.3538595939086505e-02, -1.6891587875325020e-02,
    -1.8646984919441702e-02, -1.8149697899123560e-02, -1.4875876924586697e-02,
    -8.5110608557150517e-03, 9.8848931927316319e-04,  1.3360421141947857e-02,
    2.8033301291042201e-02,  4.4158668590312596e-02,  6.0676486642862550e-02,
    7.6408062643700314e-02,  9.0166807112971648e-02,  1.0087463525509034e-01,
    1.0767073207825099e-01,  1.1000000000000000e-01,  1.0767073207825099e-01,
    1.0087463525509034e-01,  9.0166807112971648e-02,  7.6408062643700314e-02,
    6.0676486642862550e-02,  4.4158668590312596e-02,  2.8033301291042201e-02,
    1.3360421141947857e-02,  9.8848931927316319e-04,  -8.5110608557150517e-03,
    -1.4875876924586697e-02, -1.8149697899123560e-02, -1.8646984919441702e-02,
    -1.6891587875325020e-02, -1.3538595939086505e-02, -9.2899200683230643e-03,
    -4.8141100026888716e-03, -6.7950808883015244e-04, 2.6926388909554420e-03,
    5.0516421324586546e-03,  6.3172870282586493e-03,  6.5591485566565593e-03,
    5.9606317814670249e-03,  4.7744814146589041e-03,  3.2768418720744217e-03,
    1.7262548802593762e-03,  3.3245249132384837e-04,  -7.6225514924622853e-04,
    -1.4884258721727399e-03, -1.8432363328817916e-03, -1.8755700366336781e-03,
    -1.6663618496969908e-03, -1.3082618178394971e-03, -8.8813878356223768e-04,
    -4.7488111478632532e-04, -1.1360489934931548e-04, 1.7398667197948182e-04,
    3.8449091272701525e-04,  5.2507143315267811e-04,  6.0546138291252597e-04};
//initialize array for the iir a coefficients
const static double
    iir_a_coeffs[FILTER_FREQUENCY_COUNT][FILTER_FREQUENCY_COUNT] = {
        {-5.9637727070164015e+00, 1.9125339333078248e+01,
         -4.0341474540744173e+01, 6.1537466875368821e+01,
         -7.0019717951472188e+01, 6.0298814235238872e+01,
         -3.8733792862566290e+01, 1.7993533279581058e+01,
         -5.4979061224867651e+00, 9.0332828533799547e-01},
        {-4.6377947119071443e+00, 1.3502215749461572e+01,
         -2.6155952405269755e+01, 3.8589668330738348e+01,
         -4.3038990303252632e+01, 3.7812927599537133e+01,
         -2.5113598088113793e+01, 1.2703182701888094e+01,
         -4.2755083391143520e+00, 9.0332828533800291e-01},
        {-3.0591317915750960e+00, 8.6417489609637634e+00,
         -1.4278790253808875e+01, 2.1302268283304372e+01,
         -2.2193853972079314e+01, 2.0873499791105537e+01,
         -1.3709764520609468e+01, 8.1303553577932188e+00,
         -2.8201643879900726e+00, 9.0332828533800769e-01},
        {-1.4071749185996747e+00, 5.6904141470697471e+00,
         -5.7374718273676217e+00, 1.1958028362868873e+01,
         -8.5435280598354382e+00, 1.1717345583835918e+01,
         -5.5088290876998407e+00, 5.3536787286077372e+00,
         -1.2972519209655518e+00, 9.0332828533799414e-01},
        {8.2010906117760141e-01, 5.1673756579268559e+00, 3.2580350909220819e+00,
         1.0392903763919172e+01, 4.8101776408668879e+00, 1.0183724507092480e+01,
         3.1282000712126603e+00, 4.8615933365571822e+00, 7.5604535083144497e-01,
         9.0332828533799658e-01},
        {2.7080869856154512e+00, 7.8319071217995688e+00, 1.2201607990980744e+01,
         1.8651500443681620e+01, 1.8758157568004549e+01, 1.8276088095999022e+01,
         1.1715361303018897e+01, 7.3684394621253499e+00, 2.4965418284511904e+00,
         9.0332828533800436e-01},
        {4.9479835250075892e+00, 1.4691607003177602e+01, 2.9082414772101060e+01,
         4.3179839108869331e+01, 4.8440791644688879e+01, 4.2310703962394342e+01,
         2.7923434247706432e+01, 1.3822186510471010e+01, 4.5614664160654357e+00,
         9.0332828533799958e-01},
        {6.1701893352279846e+00, 2.0127225876810336e+01, 4.2974193398071684e+01,
         6.5958045321253451e+01, 7.5230437667866596e+01, 6.4630411355739852e+01,
         4.1261591079244127e+01, 1.8936128791950534e+01, 5.6881982915180291e+00,
         9.0332828533799803e-01},
        {7.4092912870072398e+00, 2.6857944460290135e+01, 6.1578787811202247e+01,
         9.8258255839887312e+01, 1.1359460153696298e+02, 9.6280452143026082e+01,
         5.9124742025776392e+01, 2.5268527576524203e+01, 6.8305064480743081e+00,
         9.0332828533799969e-01},
        {8.5743055776347692e+00, 3.4306584753117889e+01, 8.4035290411037053e+01,
         1.3928510844056814e+02, 1.6305115418161620e+02, 1.3648147221895786e+02,
         8.0686288623299745e+01, 3.2276361903872115e+01, 7.9045143816244696e+00,
         9.0332828533799636e-01}};
//initialize the array for the b coefficients 
const static double iir_b_coeffs[FILTER_FREQUENCY_COUNT][IIR_COEF_COUNT] = {
    {9.0928451882350956e-10, -0.0000000000000000e+00, -4.5464225941175478e-09,
     -0.0000000000000000e+00, 9.0928451882350956e-09, -0.0000000000000000e+00,
     -9.0928451882350956e-09, -0.0000000000000000e+00, 4.5464225941175478e-09,
     -0.0000000000000000e+00, -9.0928451882350956e-10},
    {9.0928689481426145e-10, 0.0000000000000000e+00, -4.5464344740713075e-09,
     0.0000000000000000e+00, 9.0928689481426151e-09, 0.0000000000000000e+00,
     -9.0928689481426151e-09, 0.0000000000000000e+00, 4.5464344740713075e-09,
     0.0000000000000000e+00, -9.0928689481426145e-10},
    {9.0928675031518918e-10, 0.0000000000000000e+00, -4.5464337515759461e-09,
     0.0000000000000000e+00, 9.0928675031518922e-09, 0.0000000000000000e+00,
     -9.0928675031518922e-09, 0.0000000000000000e+00, 4.5464337515759461e-09,
     0.0000000000000000e+00, -9.0928675031518918e-10},
    {9.0928691911265348e-10, 0.0000000000000000e+00, -4.5464345955632681e-09,
     0.0000000000000000e+00, 9.0928691911265362e-09, 0.0000000000000000e+00,
     -9.0928691911265362e-09, 0.0000000000000000e+00, 4.5464345955632681e-09,
     0.0000000000000000e+00, -9.0928691911265348e-10},
    {9.0928650773631441e-10, 0.0000000000000000e+00, -4.5464325386815723e-09,
     0.0000000000000000e+00, 9.0928650773631445e-09, 0.0000000000000000e+00,
     -9.0928650773631445e-09, 0.0000000000000000e+00, 4.5464325386815723e-09,
     0.0000000000000000e+00, -9.0928650773631441e-10},
    {9.0928643877976615e-10, -0.0000000000000000e+00, -4.5464321938988311e-09,
     -0.0000000000000000e+00, 9.0928643877976621e-09, -0.0000000000000000e+00,
     -9.0928643877976621e-09, -0.0000000000000000e+00, 4.5464321938988311e-09,
     -0.0000000000000000e+00, -9.0928643877976615e-10},
    {9.0928309627600387e-10, -0.0000000000000000e+00, -4.5464154813800196e-09,
     -0.0000000000000000e+00, 9.0928309627600391e-09, -0.0000000000000000e+00,
     -9.0928309627600391e-09, -0.0000000000000000e+00, 4.5464154813800196e-09,
     -0.0000000000000000e+00, -9.0928309627600387e-10},
    {9.0929565913508067e-10, 0.0000000000000000e+00, -4.5464782956754031e-09,
     0.0000000000000000e+00, 9.0929565913508061e-09, 0.0000000000000000e+00,
     -9.0929565913508061e-09, 0.0000000000000000e+00, 4.5464782956754031e-09,
     0.0000000000000000e+00, -9.0929565913508067e-10},
    {9.0926542439871930e-10, 0.0000000000000000e+00, -4.5463271219935958e-09,
     0.0000000000000000e+00, 9.0926542439871916e-09, 0.0000000000000000e+00,
     -9.0926542439871916e-09, 0.0000000000000000e+00, 4.5463271219935958e-09,
     0.0000000000000000e+00, -9.0926542439871930e-10},
    {9.0907017986989547e-10, 0.0000000000000000e+00, -4.5453508993494776e-09,
     0.0000000000000000e+00, 9.0907017986989551e-09, 0.0000000000000000e+00,
     -9.0907017986989551e-09, 0.0000000000000000e+00, 4.5453508993494776e-09,
     0.0000000000000000e+00, -9.0907017986989547e-10}};
//array declaration for the power values for each output
static double currentPowerValue[POWER_VAL_SIZE];

// initializes all attributes to zero
void initXQueue() {
  // initialize queue
  queue_init(&(xQueue), X_QUEUE_SIZE, "X QUEUE");
  // for each value in the queue push a zero
  for (uint32_t j = 0; j < X_QUEUE_SIZE; j++) {
    queue_overwritePush(&(xQueue), 0.0);
  }
}

// initializes all attributes to zero
void initYQueue() {
  // initialize queue
  queue_init(&(yQueue), Y_QUEUE_SIZE, "Y QUEUE");
  // for each value in the queue push a zero
  for (uint32_t j = 0; j < Y_QUEUE_SIZE; j++) {
    queue_overwritePush(&(yQueue), 0.0);
  }
}
// initializes all attributes to zero
void initZQueues() {
  // iterate through each queue and update values
  for (uint32_t i = 0; i < FILTER_FREQUENCY_COUNT; i++) {
    char *newName[DEFAULT_NAME_SIZE];
    sprintf(*newName, "Z Queue: %d", i);
    queue_init(&(zQueue[i]), Z_QUEUE_SIZE, *newName);
    // for each element in queue set to zero
    for (uint32_t j = 0; j < Z_QUEUE_SIZE; j++) {
      queue_overwritePush(&(zQueue[i]), 0.0);
    }
  }
}
// initializes all attributes to zero
void initOutputQueues() {
  // iterate through each queue and update values
  for (uint32_t i = 0; i < FILTER_FREQUENCY_COUNT; i++) {
    char *newName[DEFAULT_NAME_SIZE];
    sprintf(*newName, "Output Queue: %d", i);
    queue_init(&(outputQueue[i]), OUTPUT_QUEUE_SIZE, *newName);
    // for each element in queue set to zero
    for (uint32_t j = 0; j < OUTPUT_QUEUE_SIZE; j++) {
      queue_overwritePush(&(outputQueue[i]), 0.0);
    }
  }
}

// Must call this prior to using any filter functions.
void filter_init() {
  // Init queues and fill them with 0s.
  initXQueue();       // Call queue_init() on xQueue and fill it with zeros.
  initYQueue();       // Call queue_init() on yQueue and fill it with zeros.
  initZQueues();      // Call queue_init() on all of the zQueues and fill each z
                      // queue with zeros.
  initOutputQueues(); // Call queue_init() all of the outputQueues and fill each
  // outputQueue with zeros.
}

// Use this to copy an input into the input queue of the FIR-filter (xQueue).
void filter_addNewInput(double x) {
  // push a value on to the x Queue
  queue_overwritePush(&(xQueue), x);
}

// initialize the queue values in a given queue
void filter_fillQueue(queue_t *q, double fillValue) {
  // iterate over each element in the queue q
  for (queue_size_t j = 0; j < q->size; j++) {
    // after executing this function, the queue will contain 10 values
    for (uint32_t j = 0; j < q->size; j++) {
      queue_overwritePush(q, fillValue);
    }
  }
}

// Invokes the FIR-filter. Input is contents of xQueue.
// Output is returned and is also pushed on to yQueue.
double filter_firFilter() {
  double tempZ = 0.0;
  // iterating through every value of the input data and then calculating
  // transfer function
  for (uint32_t j = 0; j < FIR_COEF_COUNT; j++) {
    tempZ +=
        queue_readElementAt(&(xQueue), FIR_COEF_COUNT - 1 - j) * fir_coeffs[j];
  }
  // overwrite the oldest element with the new filtered value
  queue_overwritePush(&(yQueue), tempZ);

  return tempZ;
}

// Use this to invoke a single iir filter. Input comes from yQueue.
// Output is returned and is also pushed onto zQueue[filterNumber].
double filter_iirFilter(uint16_t filterNumber) {

  double tempb = 0.0;
  double tempa = 0.0;
  double tempz = 0.0;
  // iterating through every value of the input data and then calculating
  // transfer function
  for (uint32_t i = 0; i < Y_QUEUE_SIZE; i++) {
    tempb += queue_readElementAt(&yQueue, Y_QUEUE_SIZE - 1 - i) *
             iir_b_coeffs[filterNumber][i];
  }
  // iterate through each element in a coefficients and multtiply them by the
  // element in 7 y queue and then add to tempa
  for (uint32_t i = 0; i < FILTER_FREQUENCY_COUNT; i++) {
    tempa += queue_readElementAt(&zQueue[filterNumber], Z_QUEUE_SIZE - 1 - i) *
             iir_a_coeffs[filterNumber][i];
  }

  // output z is equal to temp b minus temp a
  tempz = tempb - tempa;

  // overwrite the oldest element with the new filtered value in both output
  // queue and zqueue
  queue_overwritePush(&outputQueue[filterNumber], tempz);
  queue_overwritePush(&zQueue[filterNumber], tempz);
  return tempz;
}

// Use this to compute the power for values contained in an outputQueue.
// If force == true, then recompute power by using all values in the
// outputQueue. This option is necessary so that you can correctly compute
// power values the first time. After that, you can incrementally compute
// power values by:
// 1. Keeping track of the power computed in a previous run, call this
// prev-power.
// 2. Keeping track of the oldest outputQueue value used in a previous run,
// call this oldest-value.
// 3. Get the newest value from the power queue, call this newest-value.
// 4. Compute new power as: prev-power - (oldest-value * oldest-value) +
// (newest-value * newest-value). Note that this function will probably need
// an array to keep track of these values for each of the 10 output queues.
// oldest value is used to calculate the oldest value in a certain filter queue
double oldest_value[FILTER_FREQUENCY_COUNT] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0, 0.0, 0.0};
//newest value is used to store and access the newest value in a filter queue
double newest_value[FILTER_FREQUENCY_COUNT] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0, 0.0, 0.0};

//compute power will use the filter number and the output queue to calculate the power comming from each of the filters 
//to determine which player frequency is recieved
double filter_computePower(uint16_t filterNumber, bool forceComputeFromScratch,
                           bool debugPrint) {
  double sum = 0;

  // iterate through each value in the filter and compute its energy by
  // squaring it and adding to the next value
  if (forceComputeFromScratch) {
    //iterate through each of the elements in the output queue for a certain filter and then sum their squares
    for (int32_t i = 0; i < queue_size(&outputQueue[filterNumber]); i++) {
      sum += queue_readElementAt(&outputQueue[filterNumber], i) *
             queue_readElementAt(&outputQueue[filterNumber], i);
    }
    currentPowerValue[filterNumber] = sum;

  }  else {
    //if not computing from scratch, use the newest and oldest value to compute power
    newest_value[filterNumber] =
        queue_readElementAt(&outputQueue[filterNumber],
                            (queue_size(&outputQueue[filterNumber]) - 1));
    sum = currentPowerValue[filterNumber] -
          (oldest_value[filterNumber] * oldest_value[filterNumber]) +
          (newest_value[filterNumber] * newest_value[filterNumber]);

    currentPowerValue[filterNumber] = sum;
  }
  //compute oldest value
  oldest_value[filterNumber] =
      queue_readElementAt(&outputQueue[filterNumber], 0);
  return currentPowerValue[filterNumber];
}

// Returns the last-computed output power value for the IIR filter
// [filterNumber].
double filter_getCurrentPowerValue(uint16_t filterNumber) {
  return filter_computePower(filterNumber, false, false);
}

// Get a copy of the current power values.
// This function copies the already computed values into a previously-declared
// array so that they can be accessed from outside the filter software by the
// detector. Remember that when you pass an array into a C function, changes
// to the array within that function are reflected in the returned array.
void filter_getCurrentPowerValues(double powerValues[]) {
  // run through each element in the current power array and update power values
  for (uint16_t i = 0; i < POWER_VAL_SIZE; i++) {
    powerValues[i] = filter_getCurrentPowerValue(i);
  }
}

//get the normalized values for each power value in the power value array
void filter_getNormalizedPowerValues(double normalizedArray[],
                                     uint16_t *indexOfMaxValue) {
  // iterate through the current power value array and find the largest value
  // divide each value by the max value and store in a normalized array
  for (int32_t i = 0; i < POWER_VAL_SIZE; i++) {
    normalizedArray[i] =
        filter_getCurrentPowerValue(i) / currentPowerValue[*indexOfMaxValue];
  }
}

/*********************************************************************************************************
********************************** Verification-assisting functions.
**************************************
********* Test functions access the internal data structures of the filter.c
*via these functions. ********
*********************** These functions are not used by the main filter
*functions. ***********************
**********************************************************************************************************/

// Returns the array of FIR coefficients.
const double *filter_getFirCoefficientArray() { return fir_coeffs; }

// Returns the number of FIR coefficients.
uint32_t filter_getFirCoefficientCount() { return FIR_COEF_COUNT; }

// Returns the array of coefficients for a particular filter number.
const double *filter_getIirACoefficientArray(uint16_t filterNumber) {
  return &iir_a_coeffs[filterNumber][0];
}

// Returns the number of A coefficients.
uint32_t filter_getIirACoefficientCount() { return IIR_A_COEF_COUNT; }

// Returns the array of coefficients for a particular filter number.
const double *filter_getIirBCoefficientArray(uint16_t filterNumber) {
  return &iir_b_coeffs[filterNumber][0];
}

// Returns the number of B coefficients.
uint32_t filter_getIirBCoefficientCount() { return IIR_COEF_COUNT; }

// Returns the size of the yQueue.
uint32_t filter_getYQueueSize() { return Y_QUEUE_SIZE; }

// Returns the decimation value.
uint16_t filter_getDecimationValue() { return DECIMATION_FACTOR; }

// Returns the address of xQueue.
queue_t *filter_getXQueue() { return &xQueue; }

// Returns the address of yQueue.
queue_t *filter_getYQueue() { return &yQueue; }

// Returns the address of zQueue for a specific filter number.
queue_t *filter_getZQueue(uint16_t filterNumber) {
  return &zQueue[filterNumber];
}

// Returns the address of the IIR output-queue for a specific filter-number.
queue_t *filter_getIirOutputQueue(uint16_t filterNumber) {
  return &outputQueue[filterNumber];
}
// void filter_runTest();
