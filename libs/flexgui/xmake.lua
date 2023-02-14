add_requires("glm")

target("flexgui")
   set_kind("static")
   add_packages("imgui", "glm")
   add_deps("implot")
   add_files("src/*.cpp")
   add_includedirs("include", {public = true})
