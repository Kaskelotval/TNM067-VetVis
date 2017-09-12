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

#ifndef IVW_MARCHINGTETRAHEDRA_H
#define IVW_MARCHINGTETRAHEDRA_H

#include <modules/tnm067lab2/tnm067lab2moduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/ports/meshport.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>

namespace inviwo {

class IVW_MODULE_TNM067LAB2_API MarchingTetrahedra : public Processor { 
public:
    struct  HashFunc
    {
        static size_t max;
        size_t operator()(std::pair<size_t,size_t> p) const{
            size_t asdf;
            asdf = p.first;
            asdf += p.second*max;
            return std::hash<size_t>{}(asdf);
        }
    };


    struct Voxel {
        vec3 pos;
        float value;
        size_t index;
    };

    struct Cell {
        Voxel voxels[8];
    };

    struct Tetrahedra {
        Voxel voxels[4];
    };


    struct MeshHelper {

        MeshHelper(std::shared_ptr<const Volume> vol);

        /**
         * Adds a vertex to the mesh. The input parameters i and j are the voxel-indices of the two
         * voxels spanning the edge on which the vertex lies. The vertex will only be added created
         * if a vertex between the same indices has not been added before. Will return the index of
         * the created vertex or the vertex that was created for this edge before. The voxel-index i
         * and j can be given in any order. 
         *
         * @param pos spatial position of the vertex
         * @param i voxel index of first voxel of the edge
         * @param j voxel index of second voxel of the edge
         */
        std::uint32_t addVertex(vec3 pos, size_t i, size_t j);
        void addTriangle(size_t i0, size_t i1, size_t i2);
        std::shared_ptr<BasicMesh> toBasicMesh();

    private:


        std::unordered_map<std::pair<size_t, size_t>, size_t , HashFunc > edgeToVertex_;
        std::vector<BasicMesh::Vertex> vertices_;
        std::shared_ptr<BasicMesh> mesh_;
        IndexBufferRAM* indexBuffer_;
    };


    MarchingTetrahedra();
    virtual ~MarchingTetrahedra() = default;
     
    virtual void process() override;

    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;
private:
    VolumeInport volume_;
    MeshOutport mesh_;

    FloatProperty isoValue_;
};

} // namespace

#endif // IVW_MARCHINGTETRAHEDRA_H

