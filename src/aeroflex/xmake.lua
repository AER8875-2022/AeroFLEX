target("aeroflex")
    set_kind("binary")

    add_deps("flexgui", "implot", "rans", "vlm")

    add_packages("imgui", "eigen", "openmp")
    add_files("src/*.cpp")
    add_includedirs("../common")
