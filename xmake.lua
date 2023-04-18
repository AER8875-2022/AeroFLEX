set_project("AeroFLEX")
set_version("0.1.0")
set_license("GPL-2.0")

add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_warnings("all")
set_languages("c++17")

if is_mode("debug") and is_plat("windows") then
    add_cxxflags("/bigobj")
end

add_requires("openmp", "eigen")
add_requires("imgui v1.89-docking", {configs = {glfw_opengl3 = true}})
add_includedirs("libs/headeronly")
includes("libs/**/xmake.lua")
includes("src/**/xmake.lua")