#pragma once

bool g_ApplicationRunning = true;

namespace FlexGUI {
	extern int Main(int argc, char** argv);
}

// #ifdef _WIN32

// #define WIN32_LEAN_AND_MEAN
// #define NOMINMAX
// #include <Windows.h>

// int APIENTRY WinMain(HINSTANCE hInst, HINSTANCE hInstPrev, PSTR cmdline, int cmdshow)
// {
// 	return FlexGUI::Main(__argc, __argv);
// }
// #else

// int main(int argc, char** argv)
// {
// 	return FlexGUI::Main(argc, argv);
// }

// #endif // _WIN32

int main(int argc, char** argv)
{
	return FlexGUI::Main(argc, argv);
}