target("implot")
   set_kind("static")
   add_packages("imgui")
   add_files("src/*.cpp")
   add_includedirs("include", {public = true})
