#include "field.h"
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/classes/surface_tool.hpp>

using namespace godot;
using namespace Eigen;

void Field::_bind_methods() {
}

Field::Field() : n(8), xt(n, n, n) {
    xt.setRandom();
    yt = xt.random();
    zt = xt.random();
    xy = xt.random();
    yz = xt.random();
    zx = xt.random();
}

void Field::_ready() {
    SurfaceTool st;
    st.begin(Mesh::PRIMITIVE_LINES);

    st.add_vertex(Vector3(0, 0, 0));
    st.add_vertex(Vector3(1, 0, 0));
    st.add_vertex(Vector3(5./12, 1./6, 0));
    st.add_vertex(Vector3(7./12, 0, 0));
    st.add_vertex(Vector3(7./12, 0, 0));
    st.add_vertex(Vector3(5./12, -1./6, 0));

    auto arrow = st.commit();
    st.clear();
    PackedColorArray colors;

    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            for (int z = 0; z < n; z++) {
                PackedColorArray color;
                color.resize(6);
                color.fill(Color(xt(x, y, z), 1., 0.));
                colors.append_array(color);
                
                st.append_from(arrow, 0, Transform3D(Basis(), Vector3(x, y, z)));
    	    }
    	}
    }

    Array arrays = st.commit_to_arrays();
    arrays[Mesh::ARRAY_COLOR] = colors;
    Ref<ArrayMesh> mesh;
    mesh.instantiate();
    mesh->add_surface_from_arrays(Mesh::PRIMITIVE_LINES, arrays);
    set_mesh(mesh);
}

void Field::_process(double delta) {
}
