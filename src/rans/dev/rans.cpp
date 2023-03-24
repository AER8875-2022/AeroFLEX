#include <rans/rans.h>

using namespace rans;

int main(int argc, char **argv) {
    tiny::config io;
    GUIHandler dummy_gui;
	bool success = io.read("../../../../examples/rans/config.ini");
	if (!success) return 1;

    Rans rans(dummy_gui);
    rans.settings.import_config_file(io);
    rans.input();
    rans.solve();

    return 0;
}
