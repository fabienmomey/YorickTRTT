# YorickTRTT
TRTT : iTeraTive Reconstruction for Tomography => code developed with Yorick for 2D/3D/4D tomographic reconstruction using inverse approach

# Quick description
TRTT is a code used for creating algebraic operators for performing Radon transform for Tomography. It is based on the modelization of the voxelized tomographic object on a basis of B-splines functions, approximated to spherically symmetric basis functions (blobs).

TRTT is implemented with Yorick [http://yorick.sourceforge.net/]. It requires several Yorick tools implemented by Éric Thiébaut for signal processing and inverse problems. These tools are available on Github: [https://github.com/emmt]. 

DO NOT HESITATE TO VISIT, "LIKE" AND "FOLLOW" THIS PAGE ;D !!!

# Requirements

Clone and install the following repositories :

- OptimPack (convex optimization algorithms): 

  		git clone https://github.com/emmt/OptimPack.git

- YTotvar (Edge-preserving regularization for Yorick): 
  		  
		git clone https://github.com/emmt/YTotVar.git

- SPL (Signal Processing Library): 
  		  
		available soon !

# To start
Get TRTT plugin:

		git clone https://github.com/fabienmomey/YorickTRTT.git

Compilation of the TRTT plugin:
	     
		$ cd YoritckTRTT/

		$ cd plugin/

		$ yorick -batch make.i

		$ make

An additional plugin is available, containing the implementation of 2D projectors by the Long & Fessler method [Long, Fessler and Balter, 2010, 3D Forward and Back-Projection for X-Ray CT Using Separable Footprints]. 

   	      	$ cd tomography-LongFessler-model/

		$ yorick -batch make.i

		$ make

		$ make install
	
To use the code, go to YoritckTRTT/ directory. You have to create en environment variable TRTT_PATH specifying the path of the YorickTRTT directory:

       	   	$ export TRTT_PATH=$(pwd) (if you are at YorickTRTT/)

Go reconstruct hot sources on cold backgrounds and planar movements ;D !!! (Tribute to PG)
