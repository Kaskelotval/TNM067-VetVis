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

#include <modules/tnm067common/tnm067commonmodule.h>
#include <modules/tnm067common/processors/gauss2dfunction.h>
#include <modules/tnm067common/processors/test2by2image.h>

namespace inviwo {

TNM067CommonModule::TNM067CommonModule(InviwoApplication* app) : InviwoModule(app, "TNM067Common") {   
    registerProcessor<Gauss2DFunction>();
    registerProcessor<Test2by2Image>();
    // Add a directory to the search path of the Shadermanager
    // ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register objects that can be shared with the rest of inviwo here:
    
    // Processors
    // registerProcessor<TNM067CommonProcessor>();
    
    // Properties
    // registerProperty<TNM067CommonProperty>();
    
    // Readers and writes
    // registerDataReader(util::make_unique<TNM067CommonReader>());
    // registerDataWriter(util::make_unique<TNM067CommonWriter>());
    
    // Data converters
    // registerRepresentationConverter(util::make_unique<TNM067CommonDisk2RAMConverter>());

    // Ports
    // registerPort<TNM067CommonOutport>("TNM067CommonOutport");
    // registerPort<TNM067CommonInport>("TNM067CommonInport");

    // PropertyWidgets
    // registerPropertyWidget<TNM067CommonPropertyWidget, TNM067CommonProperty>("Default");
    
    // Dialogs
    // registerDialog<TNM067CommonDialog>(TNM067CommonOutport);
    
    // Other varius things
    // registerCapabilities(util::make_unique<TNM067CommonCapabilities>());
    // registerSettings(util::make_unique<TNM067CommonSettings>());
    // registerMetaData(util::make_unique<TNM067CommonMetaData>());   
    // registerPortInspector("TNM067CommonOutport", "path/workspace.inv");
    // registerProcessorWidget(std::string processorClassName, std::unique_ptr<ProcessorWidget> processorWidget);
    // registerDrawer(util::make_unique_ptr<TNM067CommonDrawer>());  
}

} // namespace
