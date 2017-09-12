/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2016 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 *FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#include <inviwo/core/util/logcentral.h>
#include <modules/opengl/texture/textureutils.h>
#include <modules/tnm067lab1/processors/imageupsampler.h>
#include <modules/tnm067lab1/utils/interpolationmethods.h>
#include <inviwo/core/datastructures/image/layerram.h>
#include <inviwo/core/datastructures/image/layerramprecision.h>
#include <inviwo/core/util/imageramutils.h>

namespace inviwo {

    namespace detail{

        


        template<typename T>
        void upsample( ImageUpsampler::IntepolationMethod  method , 
            const LayerRAMPrecision<T> &inputImage, LayerRAMPrecision<T> &outputImage){
            using F = typename float_type<T>::type;
            
            size2_t inputSize = inputImage.getDimensions();
            size2_t outputSize = outputImage.getDimensions();

            const T* inPixels = inputImage.getDataTyped();
            T* outPixels = outputImage.getDataTyped();


            auto inIndex = [&inputSize](auto pos) -> size_t{
                pos = glm::clamp(pos, decltype(pos)(0) , decltype(pos)(inputSize - 1ull));
                return pos.x + pos.y * inputSize.x;
            };
            auto outIndex = [&outputSize](auto pos) -> size_t{
                pos = glm::clamp(pos, decltype(pos)(0) , decltype(pos)(outputSize - 1ull));
                return pos.x + pos.y * outputSize.x;
            };

            util::forEachPixel(outputImage, [&](ivec2 outImageCoords) {
                // outImageCoords: Exact pixel coordinates in the output image currently writing to
                // inImageCoords: Relative coordinates of outImageCoords in the input image, might be between pixels
                dvec2 inImageCoords = ImageUpsampler::convertCoordinate(outImageCoords, inputSize, outputSize);

                T finalColor(0);

                //DUMMY COLOR, remove or overwrite this bellow
                finalColor = inPixels[inIndex( glm::clamp(size2_t(outImageCoords), size2_t(0), size2_t(inputSize - size2_t(1))) )];


                // TASK 12 : Implement the following functions using TNM067::Interpolation::bilinear etc.  
                switch (method) {
                    case inviwo::ImageUpsampler::IntepolationMethod::PiecewiseConstant:
					{
						// Update finalColor
						finalColor = inPixels[inIndex(ivec2(std::floor(inImageCoords.x), std::floor(inImageCoords.y)))];
						break;
                    }
                    case inviwo::ImageUpsampler::IntepolationMethod::Bilinear:
                    {
						//						finalColor = inviwo::TNM067::Interpolation::linear(inPixels[inIndex(ivec2(floor(inImageCoords.x),floor(inImageCoords.y))]), inPixels[inIndex(ivec2(floor(outImageCoords.x), floor(outImageCoords.y))]));
						std::array<T, 4> a = {
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 1	, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y) + 1))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 1	, std::floor(inImageCoords.y) + 1))]
						};
							finalColor = inviwo::TNM067::Interpolation::bilinear(
								a,
							    inImageCoords.x - std::floor(inImageCoords.x),
								inImageCoords.y - std::floor(inImageCoords.y));
						// Update finalColor
                        break;
                    }
                    case inviwo::ImageUpsampler::IntepolationMethod::Quadratic:
                    {
                        // Update finalColor
						std::array<T, 9> a = {
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 1	, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 2	, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y) +1))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)	+1	, std::floor(inImageCoords.y) +1))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)	+2	, std::floor(inImageCoords.y)+1))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y) +2))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 1	, std::floor(inImageCoords.y)+2))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 2	, std::floor(inImageCoords.y)+2))]
						};
						finalColor = inviwo::TNM067::Interpolation::biQuadratic(
							a,
							inImageCoords.x - std::floor(inImageCoords.x),
							inImageCoords.y - std::floor(inImageCoords.y));
                        
						break;
                    }
                    case inviwo::ImageUpsampler::IntepolationMethod::Barycentric:
                    {

						std::array<T, 4> a = {
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 1	, std::floor(inImageCoords.y)))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x)		, std::floor(inImageCoords.y) + 1))],
							inPixels[inIndex(ivec2(std::floor(inImageCoords.x) + 1	, std::floor(inImageCoords.y) + 1))]
						};
						finalColor = inviwo::TNM067::Interpolation::barycentric(
							a,
							inImageCoords.x - std::floor(inImageCoords.x),
							inImageCoords.y - std::floor(inImageCoords.y));
                        break;
                    }
                    default:
                        break;
                }

                outPixels[outIndex(outImageCoords)] = finalColor;
            });
        }
    
    }


const ProcessorInfo ImageUpsampler::processorInfo_{
    "org.inviwo.imageupsampler",  // Class identifier
    "Image Upsampler",            // Display name
    "TNM067",                     // Category
    CodeState::Experimental,      // Code state
    Tags::None,                   // Tags
};
const ProcessorInfo ImageUpsampler::getProcessorInfo() const { return processorInfo_; }

ImageUpsampler::ImageUpsampler()
    : Processor()
    , inport_("inport", true)
    , outport_("outport" , true)
    , interpolationMethod_(
          "interpolationMethod", "Interpolation Method",
          {
              {"piecewiseconstant", "Piecewise Constant (Nearest Neighbor)", IntepolationMethod::PiecewiseConstant},
              {"bilinear", "Bilinear", IntepolationMethod::Bilinear},
              {"quadratic", "Quadratic", IntepolationMethod::Quadratic},
              {"barycentric", "Barycentric", IntepolationMethod::Barycentric},
          }){
    addPort(inport_);
    addPort(outport_);
    addProperty(interpolationMethod_);
}

void ImageUpsampler::process() {
    auto inputImage = inport_.getData();
    if(inputImage->getDataFormat()->getComponents()!=1){
        LogError("The ImageUpsampler processor does only support single channel images");
    }

    auto inSize = inport_.getData()->getDimensions();
    auto outDim = outport_.getDimensions();

    auto outputImage = std::make_shared<Image>(outDim,inputImage->getDataFormat());
    outputImage->getColorLayer()->setSwizzleMask(inputImage->getColorLayer()->getSwizzleMask());
    outputImage->getColorLayer()->getEditableRepresentation<LayerRAM>()->dispatch<void,dispatching::filter::Scalars>([&](auto outRep){
        auto inRep = inputImage->getColorLayer()->getRepresentation<LayerRAM>();
        detail::upsample(interpolationMethod_.get(), *(const decltype(outRep))(inRep), *outRep);
    });

    outport_.setData(outputImage);
}

dvec2 ImageUpsampler::convertCoordinate(ivec2 outImageCoords, size2_t inputSize, size2_t outputsize) {
    // TODO implement
    dvec2 c(outImageCoords);
	c = c*dvec2(double(inputSize.x) / double(outputsize.x) , double(inputSize.y) / double(outputsize.y));
    // TASK 7: Convert the outImageCoords to its coordinates in the input image 

    return c; 
}



}  // namespace inviwo
