/* tables of Altschul-Gish parameters */

/* $Name: fa21_1_1 $ - $Id: alt_parms.h,v 1.2 2007/01/19 03:38:51 wrp Exp $ */


/* first entry must be for (inf,inf) penalty */

struct alt_p {
  int gap;
  int ext;
  float Lambda;
  float K;
  float H;
};

/* BL80 1/2 bit */
struct alt_p bl80_p[] = {
  {0, 0, 0.343, 0.177, 0.66},
  {14, 2, 0.336, 0.150, 0.62},
  {12, 2, 0.328, 0.130, 0.54},
  {12, 1, 0.314, 0.096, 0.41},
  {11, 2, 0.320, 0.110, 0.51},
  {11, 1, 0.296, 0.066, 0.36},
  {10, 2, 0.311, 0.097, 0.46},
  {10, 1, 0.282, 0.052, 0.29},
  { 9, 2, 0.292, 0.069, 0.33},
  { 9, 1, 0.248, 0.026, 0.18},
  { 8, 2, 0.271, 0.050, 0.27},
  { 8, 1, 0.189, 0.0071, 0.07}
};

/* BL62 1/2 bit */
struct alt_p bl62_p[] = {
  {0, 0, 0.318, 0.13, 0.40},
  {12, 3, 0.305, 0.10, 0.38},
  {12, 2, 0.300, 0.09, 0.34},
  {12, 1, 0.275, 0.05, 0.25},
  {11, 3, 0.301, 0.09, 0.36},
  {11, 2, 0.286, 0.07, 0.29},
  {11, 1, 0.255, 0.035, 0.19},
  {10, 4, 0.293, 0.08, 0.33},
  {10, 3, 0.281, 0.06, 0.29},
  {10, 2, 0.266, 0.04, 0.24},
  {10, 1, 0.216, 0.014, 0.12},
  {9, 5, 0.286, 0.08, 0.29},
  {9, 4, 0.273, 0.06, 0.25},
  {9, 4, 0.273, 0.06, 0.25},
  {9, 2, 0.244, 0.030, 0.18},
  {9, 1, 0.176, 0.008, 0.06},
  {8, 8, 0.270, 0.06, 0.25},
  {8, 7, 0.270, 0.06, 0.25},
  {8, 6, 0.262, 0.05, 0.23},
  {8, 5, 0.262, 0.05, 0.23},
  {8, 4, 0.262, 0.05, 0.23},
  {8, 3, 0.243, 0.035, 0.18},
  {8, 2, 0.215, 0.021, 0.12},
  {7, 7, 0.247, 0.05, 0.18},
  {7, 6, 0.247, 0.05, 0.18},
  {7, 5, 0.230, 0.030, 0.15},
  {7, 4, 0.230, 0.030, 0.15},
  {7, 3, 0.208, 0.021, 0.11},
  {7, 2, 0.164, 0.009, 0.06},
  {6, 6, 0.200, 0.021, 0.10},
  {6, 5, 0.200, 0.021, 0.10},
  {6, 4, 0.179, 0.014, 0.08},
  {6, 3, 0.153, 0.010, 0.05},
  {5, 5, 0.131, 0.009, 0.04},
  {-1, -1, -1.0, -1.0, -1.0},
};

/* BL50 1/3 bit */

struct alt_p bl50_p[] = {
  {0, 0, 0.232, 0.11, 0.34},
  {16, 4, 0.222, 0.08, 0.31},
  {16, 3, 0.213, 0.06, 0.27},
  {16, 2, 0.207, 0.05, 0.24},
  {16, 1, 0.180, 0.024, 0.15},
  {15, 8, 0.222, 0.09, 0.31},
  {15, 7, 0.219, 0.08, 0.29},
  {15, 6, 0.219, 0.08, 0.29},
  {15, 5, 0.216, 0.07, 0.28},
  {15, 4, 0.216, 0.07, 0.28},
  {15, 3, 0.210, 0.06, 0.25},
  {15, 2, 0.202, 0.05, 0.22},
  {15, 1, 0.166, 0.018, 0.11},
  {14, 8, 0.218, 0.08, 0.29},
  {14, 7, 0.214, 0.07, 0.27},
  {14, 6, 0.214, 0.07, 0.27},
  {14, 5, 0.214, 0.07, 0.27},
  {14, 4, 0.205, 0.05, 0.24},
  {14, 3, 0.201, 0.05, 0.22},
  {14, 2, 0.188, 0.034, 0.17},
  {14, 1, 0.140, 0.009, 0.07},
  {13, 8, 0.211, 0.06, 0.27},
  {13, 7, 0.205, 0.05, 0.24},
  {13, 6, 0.205, 0.05, 0.24},
  {13, 5, 0.205, 0.05, 0.24},
  {13, 4, 0.202, 0.05, 0.22},
  {13, 3, 0.188, 0.034, 0.18},
  {13, 2, 0.174, 0.025, 0.13},
  {13, 1, 0.114, 0.006, 0.04},
  {12, 7, 0.205, 0.06, 0.24},
  {12, 6, 0.197, 0.05, 0.21},
  {12, 5, 0.197, 0.05, 0.21},
  {12, 4, 0.192, 0.04, 0.18},
  {12, 3, 0.178, 0.028, 0.15},
  {12, 2, 0.158, 0.019, 0.10},
  {11, 8, 0.197, 0.05, 0.21},
  {11, 7, 0.190, 0.04, 0.19},
  {11, 6, 0.190, 0.04, 0.19},
  {11, 5, 0.184, 0.04, 0.17},
  {11, 4, 0.177, 0.031, 0.15},
  {11, 3, 0.167, 0.028, 0.11},
  {11, 2, 0.130, 0.009, 0.06},
  {10, 8, 0.183, 0.04, 0.17},
  {10, 7, 0.178, 0.035, 0.16},
  {10, 6, 0.178, 0.035, 0.16},
  {10, 5, 0.168, 0.026, 0.13},
  {10, 4, 0.156, 0.020, 0.10},
  {10, 3, 0.139, 0.013, 0.07},
  {10, 2, 0.099, 0.007, 0.03},
  {9, 7, 0.164, 0.029, 0.13},
  {9, 6, 0.152, 0.021, 0.10},
  {9, 5, 0.152, 0.021, 0.10},
  {9, 4, 0.134, 0.014, 0.07},
  {9, 3, 0.107, 0.008, 0.04},
  {8, 8, 0.139, 0.017, 0.08},
  {8, 7, 0.134, 0.015, 0.07},
  {8, 6, 0.127, 0.013, 0.06},
  {8, 5, 0.117, 0.011, 0.05},
  {8, 4, 0.101, 0.009, 0.03},
  {7, 7, 0.100, 0.010, 0.04},
  {7, 6, 0.094, 0.010, 0.03},
  {-1, -1, -1.0, -1.0, -1.0},
};

struct alt_p p250_p[] = {
  {0, 0, 0.229, 0.09, 0.23},
  {16, 4, 0.217, 0.07, 0.21},
  {16, 3, 0.208, 0.05, 0.18},
  {16, 2, 0.200, 0.04, 0.16},
  {16, 1, 0.172, 0.018, 0.09},
  {15, 5, 0.215, 0.06, 0.20},
  {15, 4, 0.208, 0.05, 0.18},
  {15, 3, 0.203, 0.04, 0.16},
  {15, 2, 0.193, 0.035, 0.14},
  {15, 1, 0.154, 0.012, 0.07},
  {14, 6, 0.212, 0.06, 0.19},
  {14, 5, 0.204, 0.05, 0.17},
  {14, 4, 0.204, 0.05, 0.17},
  {14, 3, 0.194, 0.035, 0.14},
  {14, 2, 0.180, 0.025, 0.11},
  {14, 1, 0.131, 0.008, 0.04},
  {13, 6, 0.206, 0.06, 0.17},
  {13, 5, 0.196, 0.04, 0.14},
  {13, 4, 0.196, 0.04, 0.14},
  {13, 3, 0.184, 0.029, 0.12},
  {13, 2, 0.163, 0.016, 0.08},
  {13, 1, 0.110, 0.008, 0.03},
  {12, 7, 0.199, 0.05, 0.15},
  {12, 6, 0.191, 0.04, 0.13},
  {12, 5, 0.191, 0.04, 0.13},
  {12, 4, 0.181, 0.029, 0.12},
  {12, 3, 0.170, 0.022, 0.10},
  {12, 2, 0.145, 0.012, 0.06},
  {11, 7, 0.186, 0.04, 0.13},
  {11, 6, 0.180, 0.034, 0.11},
  {11, 5, 0.180, 0.034, 0.11},
  {11, 4, 0.165, 0.021, 0.09},
  {11, 3, 0.153, 0.017, 0.07},
  {11, 2, 0.122, 0.009, 0.04},
  {10, 8, 0.175, 0.031, 0.11},
  {10, 7, 0.171, 0.029, 0.10},
  {10, 6, 0.165, 0.024, 0.09},
  {10, 5, 0.158, 0.020, 0.08},
  {10, 4, 0.148, 0.017, 0.07},
  {10, 3, 0.129, 0.012, 0.05},
  {9, 7, 0.151, 0.020, 0.07},
  {9, 6, 0.146, 0.019, 0.06},
  {9, 5, 0.137, 0.015, 0.05},
  {9, 4, 0.121, 0.011, 0.04},
  {9, 3, 0.102, 0.010, 0.03},
  {8, 8, 0.123, 0.014, 0.05},
  {8, 7, 0.123, 0.014, 0.05},
  {8, 6, 0.115, 0.012, 0.04},
  {8, 5, 0.107, 0.011, 0.03},
  {7, 7, 0.090, 0.014, 0.02},
  {-1, -1, -1.0, -1.0, -1.0},
};

struct alt_p p120_p[] = {
  {0, 0, 0.342, 0.19, 0.63},
  {12, 4, 0.334, 0.14, 0.60},
  {12, 3, 0.330, 0.13, 0.57},
  {12, 2, 0.330, 0.13, 0.57},
  {12, 1, 0.219, 0.11, 0.46},
  {11, 3, 0.330, 0.13, 0.57},
  {11, 2, 0.323, 0.12, 0.51},
  {11, 1, 0.296, 0.06, 0.38},
  {10, 5, 0.323, 0.12, 0.54},
  {10, 4, 0.314, 0.09, 0.50},
  {10, 3, 0.314, 0.09, 0.50},
  {10, 2, 0.301, 0.07, 0.42},
  {10, 1, 0.273, 0.04, 0.28},
  {9, 5, 0.316, 0.11, 0.49},
  {9, 4, 0.311, 0.10, 0.45},
  {9, 3, 0.311, 0.10, 0.45},
  {9, 2, 0.284, 0.05, 0.35},
  {9, 1, 0.239, 0.023, 0.18},
  {8, 6, 0.307, 0.10, 0.43},
  {8, 5, 0.295, 0.08, 0.39},
  {8, 4, 0.295, 0.08, 0.39},
  {8, 3, 0.284, 0.06, 0.34},
  {8, 2, 0.262, 0.04, 0.26},
  {8, 1, 0.183, 0.009, 0.08},
  {7, 7, 0.286, 0.08, 0.34},
  {7, 6, 0.286, 0.08, 0.34},
  {7, 5, 0.276, 0.06, 0.31},
  {7, 4, 0.276, 0.06, 0.31},
  {7, 3, 0.255, 0.04, 0.24},
  {7, 2, 0.224, 0.023, 0.16},
  {6, 6, 0.248, 0.04, 0.23},
  {6, 5, 0.248, 0.04, 0.23},
  {6, 4, 0.234, 0.033, 0.19},
  {6, 3, 0.216, 0.025, 0.15},
  {6, 2, 0.160, 0.009, 0.06},
  {5, 5, 0.191, 0.019, 0.11},
  {5, 4, 0.173, 0.013, 0.09},
  {5, 3, 0.134, 0.006, 0.05},
  {-1, -1, -1.0, -1.0, -1.0}
};

struct alt_p bl55_p[] = {
  {0, 0, 0.224, 0.12, 0.36},
  {16, 4, 0.213, 0.08, 0.32},
  {16, 3, 0.205, 0.07, 0.28},
  {16, 2, 0.198, 0.06, 0.23},
  {16, 1, 0.164, 0.020, 0.12},
  {15, 8, 0.212, 0.09, 0.31},
  {15, 7, 0.209, 0.08, 0.30},
  {15, 6, 0.209, 0.08, 0.30},
  {15, 5, 0.205, 0.07, 0.28},
  {15, 4, 0.205, 0.07, 0.28},
  {15, 3, 0.199, 0.06, 0.25},
  {15, 2, 0.190, 0.05, 0.20},
  {15, 1, 0.146, 0.013, 0.09},
  {14, 7, 0.207, 0.08, 0.29},
  {14, 6, 0.203, 0.07, 0.27},
  {14, 5, 0.203, 0.07, 0.27},
  {14, 4, 0.195, 0.05, 0.24},
  {14, 3, 0.189, 0.04, 0.21},
  {14, 2, 0.175, 0.030, 0.16},
  {14, 1, 0.119, 0.006, 0.05},
  {13, 8, 0.201, 0.07, 0.27},
  {13, 7, 0.196, 0.06, 0.24},
  {13, 6, 0.196, 0.06, 0.24},
  {13, 5, 0.196, 0.06, 0.24},
  {13, 4, 0.191, 0.05, 0.21},
  {13, 3, 0.176, 0.032, 0.17},
  {13, 2, 0.158, 0.020, 0.12},
  {12, 8, 0.195, 0.06, 0.24},
  {12, 7, 0.188, 0.05, 0.21},
  {12, 6, 0.188, 0.05, 0.21},
  {12, 5, 0.188, 0.05, 0.21},
  {12, 4, 0.180, 0.04, 0.18},
  {12, 3, 0.165, 0.026, 0.14},
  {12, 2, 0.140, 0.014, 0.08},
  {11, 8, 0.185, 0.05, 0.20},
  {11, 7, 0.179, 0.04, 0.18},
  {11, 6, 0.179, 0.04, 0.18},
  {11, 5, 0.171, 0.033, 0.16},
  {11, 4, 0.163, 0.027, 0.13},
  {11, 3, 0.151, 0.022, 0.10},
  {11, 2, 0.110, 0.008, 0.04},
  {10, 10, 0.173, 0.04, 0.16},
  {10, 9, 0.173, 0.04, 0.16},
  {10, 8, 0.167, 0.035, 0.15},
  {10, 7, 0.167, 0.035, 0.15},
  {10, 6, 0.167, 0.035, 0.15},
  {10, 5, 0.155, 0.025, 0.12},
  {10, 4, 0.142, 0.017, 0.09},
  {10, 3, 0.121, 0.011, 0.06},
  {9, 9, 0.152, 0.026, 0.11},
  {9, 8, 0.152, 0.026, 0.11},
  {9, 7, 0.152, 0.026, 0.11},
  {9, 6, 0.137, 0.018, 0.08},
  {9, 5, 0.137, 0.018, 0.08},
  {9, 4, 0.117, 0.011, 0.05},
  {9, 3, 0.090, 0.007, 0.03},
  {8, 8, 0.125, 0.014, 0.07},
  {8, 7, 0.119, 0.013, 0.06},
  {8, 6, 0.113, 0.012, 0.05},
  {8, 5, 0.102, 0.010, 0.04},
  {8, 4, 0.085, 0.009, 0.03},
  {7, 7, 0.087, 0.010, 0.03},
  {-1, -1, -1.0, -1.0, -1.0}
};

struct alt_p nt54_p[] =
{
  {0, 0, 0.192, 0.173, 0.36},
  {16, 4, 0.192, 0.177, 0.36},
  {-1, -1, -1.0, -1.0, -1.0}
};

struct alt_p rnt54_p[] =
{
  {0, 0, 0.192, 0.173, 0.36},
  {16, 4, 0.192, 0.177, 0.36},
  {-1, -1, -1.0, -1.0, -1.0}
};

struct alt_p nt32_p[] = {
  {0,  0, 0.2712, 0.131, 0.22},
  {18, 2, 0.2620, 0.100, 0.22},
  {16, 4, 0.2600, 0.098, 0.22},
  {16, 2, 0.2540, 0.081, 0.19},
  {12, 4, 0.2340, 0.054, 0.15},
  {-1, -1, -1.0, -1.0, -1.0}
};

struct alt_p nt13_p[] = {
  {0,  0,  1.374, 0.711, 1.31},
  {4,  1,  1.36,  0.67,  1.30},
  {3,  1,  1.34,  0.58,  1.19},
  {2,  1,  1.21,  0.34,  0.77},
  {-1, -1, -1.0, -1.0, -1.0}
};

/* PAM-10 (1/10 Hartley ~ 1/3 bit scale) */

struct alt_p md10_p[] = {
  {0, 0, 0.2299, 0.309, 3.45},
  {20, 4, 0.222, 0.21, 3.1},
  {20, 2, 0.218, 0.18, 2.9},
  {18, 4, 0.220, 0.20, 2.9},
  {18, 2, 0.217, 0.18, 2.7},
  {16, 4, 0.217, 0.19, 2.8},
  {16, 2, 0.212, 0.17, 2.3},
  {14, 4, 0.212, 0.17, 2.5},
  {14, 2, 0.205, 0.15, 1.9},
  {12, 4, 0.206, 0.16, 2.1},
  {12, 2, 0.190, 0.11, 1.3},
  {-1, -1, -1.0, -1.0, -1.0}
};

/* PAM-20 (1/10 Hartley ~ 1/3 bit scale) */
struct alt_p md20_p[] = {
  {0, 0, 0.230, 0.287, 2.94},
  {20, 4, 0.221, 0.19, 2.6},
  {20, 2, 0.219, 0.18, 2.5},
  {18, 4, 0.220, 0.19, 2.5},
  {18, 2, 0.218, 0.18, 2.3},
  {16, 4, 0.218, 0.18, 2.4},
  {16, 2, 0.213, 0.17, 2.0},
  {14, 4, 0.213, 0.17, 2.1},
  {14, 2, 0.204, 0.14, 1.6},
  {12, 4, 0.207, 0.17, 1.8},
  {12, 2, 0.187, 0.10, 1.1},
  {-1, -1, -1.0, -1.0, -1.0}
};

/* PAM-40 (1/10 Hartley ~ 1/3 bit scale) */
struct alt_p md40_p[] = {
  {0,  0, 0.2293, 0.257, 2.22},
  {20, 4, 0.225, 0.22, 2.1},
  {20, 2, 0.222, 0.20, 1.9},
  {18, 4, 0.224, 0.22, 2.0},
  {18, 2, 0.220, 0.20, 1.8},
  {16, 4, 0.219, 0.19, 1.8},
  {16, 2, 0.212, 0.16, 1.5},
  {14, 4, 0.211, 0.15, 1.6},
  {14, 2, 0.199, 0.11, 1.2},
  {12, 4, 0.203, 0.14, 1.3},
  {12, 2, 0.177, 0.064, 0.7},
  {-1, -1, -1.0, -1.0, -1.0}
};
