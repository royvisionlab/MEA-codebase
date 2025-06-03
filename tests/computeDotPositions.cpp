
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <math.h> 

// getXYDotTrajectories(stimFrames,motionPerFrame,spaceConstant,numDots,screenSize,seed,correlationFrames,splitContrasts)

class MexFunction : public matlab::mex::Function {
    matlab::data::ArrayFactory factory;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        const size_t numRows = inputs[0].getDimensions()[0];
        const matlab::data::TypedArray<double> positions = std::move(inputs[0]);
        const matlab::data::TypedArray<double> xyShift = std::move(inputs[1]);
        const double spaceConstant = inputs[2][0];
        
        
        int i,j;
        double distance, weight, dx, dy;
        
        
        // Generate the matrix of output positions.
        matlab::data::TypedArray<double> outShift = factory.createArray<double>({ numRows, 2 });
//         
//         for (i = 0; i < numRows; i++) {
//             outShift[i][0] = positions[i][0];
//             outShift[i][0] = positions[i][1];
//         }
        
        for (i = 0; i < numRows; i++) {
            dx = 0;
            dy = 0;
            for (j = 0; j < numRows; j++) {
                distance = computeDistance(positions[i][0], positions[i][1], positions[j][0], positions[j][1]);
                
                weight = exp(-distance / spaceConstant);
                dx += weight * xyShift[i][0];
                dy += weight * xyShift[i][1];
            }
            outShift[i][0] = dx;
            outShift[i][1] = dy;
        }
        
        outputs[0] = outShift;
    }
    
    double computeDistance(double x1, double y1, double x2, double y2)
    {
        double x = x1 - x2; //calculating number to square in next step
        double y = y1 - y2;
        double dist;

        dist = std::pow(x, 2) + std::pow(y, 2);       //calculating Euclidean distance
        dist = std::sqrt(dist);                  

        return dist;
    }
};



