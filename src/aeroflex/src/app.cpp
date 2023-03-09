#include "Application.h"
#include "EntryPoint.h"
#include "imgui.h"
#include "implot.h"
#include "tinyconfig.hpp"

#include "common_aeroflex.hpp"

#include <rans/rans.h>
#include <vlm/vlm.hpp>

#include <future>
#include <thread>
#include <optional>
#include <unordered_map>

template <class T>
bool is_future_done(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

const std::string bool_to_string(const bool b) {
    return b ? "true" : "false";
}

struct Settings {
	rans::Settings rans;
};

class Aero {
	public:
		void solve_async();
		void solve_await();
		void config_open_async(const std::string &conf_path);
		void config_open_await();
		void config_save_async(const std::string &conf_path);
		void config_save_await();

		// Modules
		rans::Rans &rans;
		vlm::VLM &vlm;

		Settings settings;

		std::future<void> future_solve;
		std::future<std::optional<Settings>> future_config_open;
		std::future<bool> future_config_save;
		std::string outfile;

		// Signal routing
		bool signal_status_busy = false;
		bool signal_status_ready = true;
		bool config_file_set = false;
		GUIHandler &gui;
		Aero(rans::Rans &rans, vlm::VLM &vlm, GUIHandler &gui);
};

struct ExampleAppLog
{
    ImGuiTextBuffer     Buf;
    ImGuiTextFilter     Filter;
    ImVector<int>       LineOffsets; // Index to lines offset. We maintain this with AddLog() calls.
    bool                AutoScroll;  // Keep scrolling if already at the bottom.

    ExampleAppLog()
    {
        AutoScroll = true;
        Clear();
    }

    void    Clear()
    {
        Buf.clear();
        LineOffsets.clear();
        LineOffsets.push_back(0);
    }

    void    AddLog(const char* fmt, ...) IM_FMTARGS(2)
    {
        int old_size = Buf.size();
        va_list args;
        va_start(args, fmt);
        Buf.appendfv(fmt, args);
        va_end(args);
        for (int new_size = Buf.size(); old_size < new_size; old_size++)
            if (Buf[old_size] == '\n')
                LineOffsets.push_back(old_size + 1);
    }

    void    Draw(const char* title, bool* p_open = NULL)
    {
        if (!ImGui::Begin(title, p_open))
        {
            ImGui::End();
            return;
        }

        // Options menu
        if (ImGui::BeginPopup("Options"))
        {
            ImGui::Checkbox("Auto-scroll", &AutoScroll);
            ImGui::EndPopup();
        }

        // Main window
        if (ImGui::Button("Options"))
            ImGui::OpenPopup("Options");
        ImGui::SameLine();
        bool clear = ImGui::Button("Clear");
        ImGui::SameLine();
        bool copy = ImGui::Button("Copy");
        ImGui::SameLine();
        Filter.Draw("Filter", -100.0f);

        ImGui::Separator();

        if (ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar))
        {
            if (clear)
                Clear();
            if (copy)
                ImGui::LogToClipboard();

            ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
            const char* buf = Buf.begin();
            const char* buf_end = Buf.end();
            if (Filter.IsActive())
            {
                // In this example we don't use the clipper when Filter is enabled.
                // This is because we don't have random access to the result of our filter.
                // A real application processing logs with ten of thousands of entries may want to store the result of
                // search/filter.. especially if the filtering function is not trivial (e.g. reg-exp).
                for (int line_no = 0; line_no < LineOffsets.Size; line_no++)
                {
                    const char* line_start = buf + LineOffsets[line_no];
                    const char* line_end = (line_no + 1 < LineOffsets.Size) ? (buf + LineOffsets[line_no + 1] - 1) : buf_end;
                    if (Filter.PassFilter(line_start, line_end))
                        ImGui::TextUnformatted(line_start, line_end);
                }
            }
            else
            {
                // The simplest and easy way to display the entire buffer:
                //   ImGui::TextUnformatted(buf_begin, buf_end);
                // And it'll just work. TextUnformatted() has specialization for large blob of text and will fast-forward
                // to skip non-visible lines. Here we instead demonstrate using the clipper to only process lines that are
                // within the visible area.
                // If you have tens of thousands of items and their processing cost is non-negligible, coarse clipping them
                // on your side is recommended. Using ImGuiListClipper requires
                // - A) random access into your data
                // - B) items all being the  same height,
                // both of which we can handle since we have an array pointing to the beginning of each line of text.
                // When using the filter (in the block of code above) we don't have random access into the data to display
                // anymore, which is why we don't use the clipper. Storing or skimming through the search result would make
                // it possible (and would be recommended if you want to search through tens of thousands of entries).
                ImGuiListClipper clipper;
                clipper.Begin(LineOffsets.Size);
                while (clipper.Step())
                {
                    for (int line_no = clipper.DisplayStart; line_no < clipper.DisplayEnd; line_no++)
                    {
                        const char* line_start = buf + LineOffsets[line_no];
                        const char* line_end = (line_no + 1 < LineOffsets.Size) ? (buf + LineOffsets[line_no + 1] - 1) : buf_end;
                        ImGui::TextUnformatted(line_start, line_end);
                    }
                }
                clipper.End();
            }
            ImGui::PopStyleVar();

            // Keep up at the bottom of the scroll region if we were already at the bottom at the beginning of the frame.
            // Using a scrollbar or mouse-wheel will take away from the bottom edge.
            if (AutoScroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
                ImGui::SetScrollHereY(1.0f);
        }
        ImGui::EndChild();
        ImGui::End();
    }
};

void solve(rans::Rans &rans) {
	rans.input();
	rans.solve();
}

std::optional<Settings> config_open(const std::string &conf_path) {
	Settings settings;

	tiny::config io;
	bool success = io.read(conf_path);
	if (!success) return std::nullopt;

	if (io.how_many("rans-bc") != 2) throw std::runtime_error("[RANS] Invalid number of boundary conditions");

	settings.rans.g.gamma = io.get<double>("rans-gas", "gamma");
	settings.rans.g.R = io.get<double>("rans-gas", "R");

	for (int i = 0; i < io.how_many("rans-bc"); i++) {
		std::string type = io.get_i<std::string>("rans-bc", "type", i);
		std::string name = io.get_i<std::string>("rans-bc", "name", i);
		settings.rans.bcs[name];
		settings.rans.bcs[name].bc_type = type;
		if (type == "farfield") {
			settings.rans.bcs[name].vars_far.T = io.get_i<double>("rans-bc", "T", i);
			settings.rans.bcs[name].vars_far.mach = io.get_i<double>("rans-bc", "mach", i);
			settings.rans.bcs[name].vars_far.angle = io.get_i<double>("rans-bc", "angle", i);
			settings.rans.bcs[name].vars_far.p = io.get_i<double>("rans-bc", "p", i);
		}
	}

	settings.rans.set_solver_type(io.get<std::string>("rans-solver", "type"));
	settings.rans.second_order = io.get<bool>("rans-solver", "second_order");
	settings.rans.relaxation = io.get<double>("rans-solver", "relaxation");

	for (int i = 0; i < io.how_many("rans-mesh"); i++) {
		settings.rans.meshes.push_back(io.get_i<std::string>("rans-mesh", "file", i));
	}
	return settings;
}

bool config_save(const std::string &conf_path, Settings &settings) {
	tiny::config io;
	io.sections = {
		"rans-gas",
		"rans-bc",
		"rans-solver",
		"rans-mesh",
	};

	io.config["rans-gas"]["gamma"] = std::to_string(settings.rans.g.gamma);
	io.config["rans-gas"]["R"] = std::to_string(settings.rans.g.R);

	io.config_vec["rans-bc"] = {};
	for (auto &[name, bc] : settings.rans.bcs) {
		io.config_vec["rans-bc"].push_back({
			{"type", bc.bc_type},
			{"name", name},
			{"T", std::to_string(bc.vars_far.T)},
			{"mach", std::to_string(bc.vars_far.mach)},
			{"angle", std::to_string(bc.vars_far.angle)},
			{"p", std::to_string(bc.vars_far.p)},
		});
	};

	io.config["rans-solver"]["type"] = settings.rans.solver_type();
	io.config["rans-solver"]["second_order"] = bool_to_string(settings.rans.second_order);
	io.config["rans-solver"]["relaxation"] = bool_to_string(settings.rans.relaxation);

	io.config_vec["rans-mesh"] = {};
	for (auto &mesh : settings.rans.meshes) {
		io.config_vec["rans-mesh"].push_back({
			{"file", mesh},
		});
	};

	return io.write(conf_path);
}

Aero::Aero(rans::Rans &rans, vlm::VLM &vlm, GUIHandler &gui) : rans(rans), vlm(vlm), gui(gui) {
	settings.rans.bcs["farfield"];
	settings.rans.bcs["farfield"].bc_type = "farfield";
	settings.rans.bcs["wall"];
	settings.rans.bcs["wall"].bc_type = "slip-wall";

	// TEMPORARY !!!!
	settings.rans.meshes.push_back("../../../../examples/rans/naca0012q_coarse.msh");
	settings.rans.meshes.push_back("../../../../examples/rans/naca0012q_mid.msh");
	settings.rans.meshes.push_back("../../../../examples/rans/naca0012q_fine.msh");
}

void Aero::solve_async() {
	signal_status_ready = false;
	signal_status_busy = true;
	rans.settings = settings.rans;
	future_solve = std::async(std::launch::async, [this](){return solve(this->rans);});
}

void Aero::solve_await() {
	future_solve.get();
	signal_status_ready = true;
	signal_status_busy = false;
	gui.signal.stop = false;
}

void Aero::config_open_async(const std::string &conf_path) {
	signal_status_ready = false;
	signal_status_busy = true;
	outfile = conf_path;
	future_config_open = std::async(std::launch::async,
	[this](const std::string &path){
		std::optional<Settings> new_settings_op{};
		try {
			new_settings_op = config_open(path);
		} catch (std::exception &e) {
			// Put this in queue
			std::cout << e.what() << std::endl;
		}
		return new_settings_op;
	}, conf_path);
}

void Aero::config_open_await() {
	auto settings_op = future_config_open.get();
	if (settings_op.has_value()) {
		config_file_set = true;
		settings = settings_op.value();
	}
	signal_status_ready = true;
	signal_status_busy = false;
}

void Aero::config_save_async(const std::string &conf_path) {
	signal_status_ready = false;
	signal_status_busy = true;
	future_config_save = std::async(std::launch::async,
	[this](const std::string &path){
		return config_save(path, this->settings);
	}, conf_path);
}

void Aero::config_save_await() {
	bool success = future_config_save.get();
	if (success) {
		gui.msg.push("Config save success\n");
	} else {
		gui.msg.push("Config save failed\n");
	}
	signal_status_ready = true;
	signal_status_busy = false;
}

struct RansLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	RansLayer(Aero &aero) : aero(aero) {};
};

struct ButtonLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	ButtonLayer(Aero &aero) : aero(aero) {};
};

struct GraphLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	GraphLayer(Aero &aero) : aero(aero) {};
};

struct ConsoleLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	ConsoleLayer(Aero &aero) : aero(aero) {};
};

void RansLayer::OnUIRender() {
	ImGui::Begin("RANS");

	ImGui::Separator();
	ImGui::Text("Gas");
	ImGui::InputDouble("gamma", &aero.settings.rans.g.gamma, 0.01f, 1.0f, "%.4f");
	ImGui::InputDouble("R", &aero.settings.rans.g.R, 0.01f, 1.0f, "%.4f");

	ImGui::Separator();
	ImGui::Text("Farfield");
	ImGui::InputDouble("Mach", &aero.settings.rans.bcs["farfield"].vars_far.mach, 0.01f, 1.0f, "%.4f");
	ImGui::InputDouble("AoA", &aero.settings.rans.bcs["farfield"].vars_far.angle, 0.01f, 1.0f, "%.4f");
	ImGui::InputDouble("Temperature", &aero.settings.rans.bcs["farfield"].vars_far.T, 0.01f, 1.0f, "%.4f");
	ImGui::InputDouble("Pressure", &aero.settings.rans.bcs["farfield"].vars_far.p, 0.01f, 1.0f, "%.4f");

	ImGui::Separator();
	ImGui::Text("Solver");

	ImGui::RadioButton("Explicit", &aero.settings.rans.type, 0); ImGui::SameLine();
	ImGui::RadioButton("Implicit", &aero.settings.rans.type, 1);

	ImGui::Checkbox("Second Order", &aero.settings.rans.second_order);
	ImGui::InputDouble("Relaxation", &aero.settings.rans.relaxation, 0.1f, 1.0f, "%.2f");

	ImGui::End();
}

void ButtonLayer::OnUIRender() {
	{
		ImGui::Begin("Buttons");
		ImGui::Text("Simulation");
		if (ImGui::Button("Start Simulation", ImVec2(-1.0f, 0.0f)) && !aero.signal_status_busy && aero.signal_status_ready) {
			aero.solve_async();
		};

		if (ImGui::Button("Stop Simulation", ImVec2(-1.0f, 0.0f)) && aero.signal_status_busy) {
			aero.gui.signal.stop = true;
		};

		if (ImGui::Button("Pause", ImVec2(-1.0f, 0.0f)) && aero.signal_status_busy) {
			aero.gui.signal.pause = true;
		};

		if (ImGui::Button("Resume", ImVec2(-1.0f, 0.0f)) && aero.signal_status_busy) {
			aero.gui.signal.pause = false;
		};

		ImGui::Separator();
		ImGui::Text("Config");
		if (ImGui::Button("Open", ImVec2(-1.0f, 0.0f)) && !aero.signal_status_busy && aero.signal_status_ready) {
			// TODO: Make this a file dialog !!
			aero.gui.msg.push("Starting parsing\n");
			aero.config_open_async("../../../../examples/conf.ini");
		};

		if (ImGui::Button("Save", ImVec2(-1.0f, 0.0f)) && aero.config_file_set) {
			aero.config_save_async("test.ini");
		};

		if (!aero.signal_status_ready && is_future_done(aero.future_solve)) {
			aero.solve_await();
		};

		if (!aero.signal_status_ready && is_future_done(aero.future_config_open)) {
			aero.config_open_await();
			// aero.gui.msg.push("Done parsing");
		};

		if (!aero.signal_status_ready && is_future_done(aero.future_config_save)) {
			aero.config_save_await();
		};

		ImGui::End();

		// ImGui::ShowDemoWindow();
		// ImPlot::ShowDemoWindow();
	};
};

void GraphLayer::OnUIRender() {
	{
		ImGui::Begin("Graphs");
		static ImPlotAxisFlags xflags = ImPlotAxisFlags_None;
		static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_RangeFit;
		const double xticks = 1;

		if (ImPlot::BeginPlot("Convergence", ImVec2(-1,400))) {
			ImPlot::SetupAxes("Iterations","Residual",xflags,yflags);
			ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, 0, aero.rans.iters);
			ImPlot::SetupAxisZoomConstraints(ImAxis_X1, 11.0, aero.rans.iters);
			ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

			if (!aero.signal_status_ready) {
				ImPlot::SetupAxisLimits(ImAxis_X1, 0, aero.rans.iters, ImPlotCond_Always);
			}
			ImPlot::PlotLine("L2 residual", aero.rans.residuals.data(), aero.rans.iters);
			ImPlot::EndPlot();
		}

		ImGui::End();
	}
};

void ConsoleLayer::OnUIRender() {
    static ExampleAppLog log;

    ImGui::Begin("Console", NULL);
    std::optional<std::string> txt = aero.gui.msg.pop();

    if (txt.has_value()) {
        log.AddLog(txt.value().c_str());
		std::cout << "should print" << std::endl;
    }
    ImGui::End();
    log.Draw("Console", NULL);

    // For the demo: add a debug button _BEFORE_ the normal log window contents
    // We take advantage of a rarely used feature: multiple calls to Begin()/End() are appending to the _same_ window.
    // Most of the contents of the window will be added by the log.Draw() call.
    // ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiCond_FirstUseEver);
    // ImGui::Begin("Example: Log", p_open);
    // if (ImGui::SmallButton("[Debug] Add 5 entries"))
    // {
    //     static int counter = 0;
    //     const char* categories[3] = { "info", "warn", "error" };
    //     const char* words[] = { "Bumfuzzled", "Cattywampus", "Snickersnee", "Abibliophobia", "Absquatulate", "Nincompoop", "Pauciloquent" };
    //     for (int n = 0; n < 5; n++)
    //     {
    //         const char* category = categories[counter % IM_ARRAYSIZE(categories)];
    //         const char* word = words[counter % IM_ARRAYSIZE(words)];
    //         log.AddLog("[%05d] [%s] Hello, current time is %.1f, here's a word: '%s'\n",
    //             ImGui::GetFrameCount(), category, ImGui::GetTime(), word);
    //         counter++;
    //     }
    // }
    // ImGui::End();

    // // Actually call in the regular Log helper (which will Begin() into the same window as we just did)
    // log.Draw("Example: Log", p_open);
}

FlexGUI::Application* CreateApplication(int argc, char** argv, Aero& aero)
{
	FlexGUI::ApplicationSpecification spec;
	spec.Name = "AEROFLEX";
	spec.Width = 1600;
	spec.Height = 900;

	FlexGUI::Application* app = new FlexGUI::Application(spec);
	app->PushLayer(std::make_shared<ButtonLayer>(aero));
	app->PushLayer(std::make_shared<GraphLayer>(aero));
	app->PushLayer(std::make_shared<RansLayer>(aero));
	app->PushLayer(std::make_shared<ConsoleLayer>(aero));

	app->SetMenubarCallback([app]()
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Exit"))
			{
				app->Close();
			}
			ImGui::EndMenu();
		}
	});
	return app;
}

namespace FlexGUI {
	int Main(int argc, char** argv)
	{
		GUIHandler gui;

		// Initialize modules with signal routing
		rans::Rans rans(gui);
		vlm::VLM vlm(gui);

		// Initialize main application with the modules
		Aero aero(rans, vlm, gui);

		while (g_ApplicationRunning)
		{
			FlexGUI::Application* app = CreateApplication(argc, argv, aero);
			app->Run();
			delete app;
		}
		return 0;
	}
}
