
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <random>

// getXYDotTrajectories(stimFrames,motionPerFrame,spaceConstant,numDots,screenSize,seed,correlationFrames,splitContrasts)

class MexFunction : public matlab::mex::Function {
    matlab::data::ArrayFactory factory;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
//         const matlab::data::TypedArray<double> Y = std::move(inputs[1]);
        const size_t stimFrames = inputs[0][0];
        const double motionPerFrame = inputs[1][0];
        const double spaceConstant = inputs[2][0];
        const size_t numDots = inputs[3][0];
        const matlab::data::TypedArray<double> screenSize = std::move(inputs[4]);
//         const int seed = 1; //(int)inputs[5][0];
//         const int * seed = const_cast<int *>(inputs[5][0]);
        std::mt19937::result_type seed = static_cast<std::mt19937::result_type>( inputs[5][0] );
        const double minRadius = inputs[6][0];
        
        // Generate the matrix of output positions.
        matlab::data::TypedArray<double> out = factory.createArray<double>({ stimFrames, numDots, 2 });
        
        int frame;
        
        for (frame = 0; frame < stimFrames; frame++) {
        }
        
        // Loop through 
        
        
        // Function implementation
//         std::random_device rd{};
        std::mt19937 gen{seed};
        
        std::normal_distribution<> d{0,1};
        
        
        out[0][0][0] = d(gen);
        
        outputs[0] = std::move(out);
    }
    
    void generateDotPositions() {
    }
};