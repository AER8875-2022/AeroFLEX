set_project("RANS")
set_version("0.1.0")

add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_warnings("all")
set_languages("c++17")

add_requires("openmp", "eigen")
add_requires("imgui v1.89-docking", {configs = {glfw_opengl3 = true}})
includes("**/xmake.lua")