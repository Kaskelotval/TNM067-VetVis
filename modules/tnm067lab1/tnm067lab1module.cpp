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

#include <modules/tnm067lab1/tnm067lab1module.h>
#include <modules/tnm067lab1/processors/imagetoheightfield.h>
#include <modules/tnm067lab1/processors/imageupsampler.h>
#include <modules/tnm067lab1/processors/imagemappingcpu.h>

namespace inviwo {

TNM067Lab1Module::TNM067Lab1Module(InviwoApplication* app) : InviwoModule(app, "TNM067Lab1") {   
    registerProcessor<ImageToHeightfield>();
    registerProcessor<ImageUpsampler>();
    registerProcessor<ImageMappingCPU>();
    // Add a directory to the search path of the Shadermanager
    // ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register objects that can be shared with the rest of inviwo here:
    
    // Processors
    // registerProcessor<TNM067Lab1Processor>();
    
    // Properties
    // registerProperty<TNM067Lab1Property>();
    
    // Readers and writes
    // registerDataReader(util::make_unique<TNM067Lab1Reader>());
    // registerDataWriter(util::make_unique<TNM067Lab1Writer>());
    
    // Data converters
    // registerRepresentationConverter(util::make_unique<TNM067Lab1Disk2RAMConverter>());

    // Ports
    // registerPort<TNM067Lab1Outport>("TNM067Lab1Outport");
    // registerPort<TNM067Lab1Inport>("TNM067Lab1Inport");

    // PropertyWidgets
    // registerPropertyWidget<TNM067Lab1PropertyWidget, TNM067Lab1Property>("Default");
    
    // Dialogs
    // registerDialog<TNM067Lab1Dialog>(TNM067Lab1Outport);
    
    // Other varius things
    // registerCapabilities(util::make_unique<TNM067Lab1Capabilities>());
    // registerSettings(util::make_unique<TNM067Lab1Settings>());
    // registerMetaData(util::make_unique<TNM067Lab1MetaData>());   
    // registerPortInspector("TNM067Lab1Outport", "path/workspace.inv");
    // registerProcessorWidget(std::string processorClassName, std::unique_ptr<ProcessorWidget> processorWidget);
    // registerDrawer(util::make_unique_ptr<TNM067Lab1Drawer>());  
}

} // namespace
