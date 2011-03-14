#ifndef MMPBSA_GLOBALS_H
#define MMPBSA_GLOBALS_H


#define MMPBSA_QUEUE_TITLE "grid_queue"//Main queue XML tag
#define MMPBSA_PI 3.14159265358979323846
#define MMPBSA_HALF_PI 1.570796326794896558
#define MMPBSA_2_PI 6.283185307179586232
#define MMPBSA_4_PI 12.566370612
#define MMPBSA_DEG_TO_RAD 0.0174532925199432957692
#define MMPBSA_FORMAT std::setprecision(5)//Precision used with importing data, e.g. Sander input files.

//Electrostatic constants
#define DEFAULT_SCNB 2.0
#define DEFAULT_SCEE 1.2
#define DEFAULT_DIELC 1.0
#define CR_CHAR 0xd

//Queue file information
#define MMPBSA_XML_TITLE "mmpbsa_energy"//Main energy data XML tag
#define SANDER_MDOUT_TYPE "mdout"
#define SANDER_PRMTOP_TYPE "prmtop"
#define MMPBSA_TOPOLOGY_TYPE "top"
#define MMPBSA_OUT_TYPE "mmpbsa_out"
#define RADII_TYPE "radii"
#define SANDER_INPCRD_TYPE "inpcrd"
#define MMPBSA_TRAJECTORY_TYPE "traj"
#define CHECKPOINT_FILE_TYPE "checkpoint"

//These tags are used to indicate the current state of the mdmmpbsa program
#define MMPBSASTATE_TAG "mmpbsa_state"
#define MDSTATE_TAG "moldyn_state"

typedef float mead_data_t;///<MEAD uses single precision. If that every changes, this provides one location to keep up with MEAD's precision.
typedef double mmpbsa_t;///<Default floating point data type. This can be used to control single versus double precision.

#endif


