/* 
 * File:   mmpbsa_graphics.h
 * Author: dcoss
 *
 * Created on August 4, 2010, 12:53 PM
 */

#ifndef MMPBSA_GRAPHICS_H
#define	MMPBSA_GRAPHICS_H

/**
 * Shared Memory is used to communicate between the science application and the
 * graphics application. MMPBSA_SHMEM stores information that is used by the
 * graphics application based on the progress of the science application.
 * 
 */
struct MMPBSA_SHMEM {
    float* crds;
    size_t natoms;
    BOINC_STATUS status;
    double update_time;
    double fraction_done;
    double cpu_time;
    

    //main program decrements this once per second.
    //graphics app sets this to 5, repeated.
    //if countdown != 0, SHMEM is updated.
    int countdown;

};

#endif	/* MMPBSA_GRAPHICS_H */

