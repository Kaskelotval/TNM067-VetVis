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
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#include <modules/tnm067common/processors/gauss2dfunction.h>
#include <inviwo/core/datastructures/image/image.h>
#include <inviwo/core/datastructures/image/layerram.h>
#include <inviwo/core/util/indexmapper.h>

namespace inviwo {

const ProcessorInfo Gauss2DFunction::processorInfo_{
    "org.inviwo.gauss2dfunction",     // Class identifier
    "Gauss2DFunction",                // Display name
    "TNM067",                         // Category
    CodeState::Experimental,          // Code state
    Tags::None,                       // Tags
};
const ProcessorInfo Gauss2DFunction::getProcessorInfo() const {
    return processorInfo_;
}

Gauss2DFunction::Gauss2DFunction()
    : Processor()
    , imageOutport_("outport", DataFloat32::get(), false)
    , size_("size", "Size", 32, 1, 256)
    , gauss2d_("gauss2d", "Gaussian") {
    
    addPort(imageOutport_);
    gauss2d_.height_ = 1.0f;
    gauss2d_.sigma_ = 0.2f;
    gauss2d_.center_ = dvec2(0.5f, 0.5f);
    addProperty(size_);
    addProperty(gauss2d_);
}
    
void Gauss2DFunction::process() {
    std::unique_ptr<Image> upScaleImage = util::make_unique<Image>(size2_t(size_.get()), DataFloat32::get());
    upScaleImage->getColorLayer()->setSwizzleMask(swizzlemasks::luminance);
    gaussImage(*static_cast<LayerRAMPrecision<float>*>(  upScaleImage->getColorLayer()->getEditableRepresentation<LayerRAM>()));
    imageOutport_.setData(upScaleImage.release());
}

void Gauss2DFunction::gaussImage(LayerRAMPrecision<float>& img) {
    auto dim = img.getDimensions();
    auto data = img.getDataTyped();
    util::IndexMapper2D index(dim);

    double height = 0.0f;
    for (size_t i = 0; i < dim.x; i++) {
        for (size_t j = 0; j < dim.y; j++) {

            double x = static_cast<double>(i) / (dim.x-1);
            double y = static_cast<double>(j) / (dim.y-1);

            //height from gauss field
            height = gauss2d_.evaluate(dvec2(x, y));
            data[index(i,j)] = height;
        }
    }
}

} // namespace

