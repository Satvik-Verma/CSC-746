//
// (C) 2021, E. Wes Bethel
// mpi_2dmesh.hpp
//
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>

using namespace std;

static char default_input_fname[] = "../data/zebra-gray-int8-4x";
static int default_data_dims[2] = {7112, 5146};
static char default_output_fname[] = "../data/processed-raw-int8-4x-cpu.dat";

typedef enum {
   ROW_DECOMP = 1,
   COLUMN_DECOMP = 2,
   TILE_DECOMP = 3
} DecompositionEnum;

typedef enum {
   NO_ACTION = 0,
   MESH_PROCESSING = 1,
   MESH_LABELING_ONLY = 2
} ActionEnum;

class AppState
{
   private:

   public:
   int myrank, nranks;
   int decomp;
   int action;
   int debug;
   char input_filename[256]; // pls update to use std::string
   char output_filename[256];  // pls update to use std::string

   int global_mesh_size[2]; // assumes 2D mesh

   // used by rank 0 to hold the float-converted input data
   vector<float> input_data_floats;

   // used to hold the results of the computation
   vector<float> output_data_floats;


   AppState(void) {
      myrank = 0;
      nranks = 1;
      decomp = ROW_DECOMP;
      // global_mesh_size[0] = global_mesh_size[1] = -1;
      global_mesh_size[0] = default_data_dims[0];
      global_mesh_size[1] = default_data_dims[1];
      debug = 0;
      action = NO_ACTION;
      strcpy(input_filename, default_input_fname);
      strcpy(output_filename, default_output_fname);

      input_data_floats.resize(0);
      output_data_floats.resize(0);
   } // end AppState constructor

   void print(void)
   {
      printf("AppState() = \n");
      printf("\tmyrank = %d \n\tnranks=%d \n\tdecomp=%d \n", myrank, nranks, decomp);
      printf("\taction = %d \n", action);
      printf("\tdims = [%d, %d] \n", global_mesh_size[0], global_mesh_size[1]);
   }

};  // class AppState

class Tile2D
{
public:
    int xloc, yloc;       // Location of the tile in the global grid (top-left corner)
    int width, height;    // Dimensions of the tile's base grid
    int ghost_xmin, ghost_xmax, ghost_ymin, ghost_ymax; // Sizes of halo regions on each side

    int tileRank;         // Rank ID owner of this tile

    vector<float> inputBuffer;  // Includes halo regions
    vector<float> outputBuffer; // Excludes halo regions

    // Constructor to initialize tile properties
    Tile2D(int tx, int ty, int xsize, int ysize, int rank)
        : xloc(tx), yloc(ty), width(xsize), height(ysize), tileRank(rank)
    {
        ghost_xmin = ghost_xmax = ghost_ymin = ghost_ymax = 0;
        inputBuffer.resize(0);  // Start with empty buffers
        outputBuffer.resize(0);
    }

    // Function to compute buffer sizes, including halos
    void allocateBuffers() {
        int paddedWidth = width + ghost_xmin + ghost_xmax;
        int paddedHeight = height + ghost_ymin + ghost_ymax;
        inputBuffer.resize(paddedWidth * paddedHeight);  // With halos
        outputBuffer.resize(width * height);            // Without halos
    }

    // Print the tile's metadata and buffers for debugging
    void print(int row, int col) {
        printf(" Tile at [%d, %d], \tx/yloc: (%d, %d),\tbase grid size [%d,%d],\trank=%d,\tgxmin/gxmax/gymin/gymax=[%d,%d,%d,%d],\tinputBuffer.size()=%d, outputBuffer.size()=%d \n",
               row, col, xloc, yloc, width, height, tileRank,
               ghost_xmin, ghost_xmax, ghost_ymin, ghost_ymax,
               inputBuffer.size(), outputBuffer.size());
    }
}; // class Tile2D

// eof
