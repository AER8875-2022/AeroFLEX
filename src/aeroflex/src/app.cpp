#include "Application.h"
#include "EntryPoint.h"
#include "imgui.h"
#include "implot.h"

#include "common_aeroflex.hpp"
#include "tinyconfig.hpp"

#include <rans/rans.h>
#include <vlm/vlm.hpp>

// Temporary
#include <rans/parser.h>

#include <future>
#include <thread>

template <class T>
bool is_future_done(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
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
		tiny::config &conf;
		rans::Rans &rans;
		vlm::VLM &vlm;

		Settings settings;

		std::future<void> future_solve;
		std::future<Settings> future_config_open;
		std::future<void> future_config_save;

		// Signal routing
		bool signal_status_busy = false;
		bool signal_status_ready = true;
		GUIHandler &gui;
		Aero(tiny::config &conf, rans::Rans &rans, vlm::VLM &vlm, GUIHandler &gui);
};

void solve(rans::Rans &rans) {
	rans.input();
	rans.solve();
}

Settings config_open(tiny::config &io, const std::string &conf_path) {
	Settings settings;
	io.read(conf_path);
	return settings;
}

Aero::Aero(tiny::config &conf, rans::Rans &rans, vlm::VLM &vlm, GUIHandler &gui) : conf(conf), rans(rans), vlm(vlm), gui(gui) {
	settings.rans.bcs["farfield"];
	settings.rans.bcs["farfield"].bc_type = "farfield";
	settings.rans.bcs["slip-wall"];
	settings.rans.bcs["slip-wall"].bc_type = "wall";

	// TEMPORARY !!!!
	settings.rans.meshes.push_back("../../../../examples/rans/naca0012q_coarse.msh");
	settings.rans.meshes.push_back("../../../../examples/rans/naca0012q_mid.msh");
	settings.rans.meshes.push_back("../../../../examples/rans/naca0012q_fine.msh");
}

void Aero::solve_async() {
	signal_status_ready = false;
	signal_status_busy = true;
	// rans.settings = settings.rans;
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
	future_config_open = std::async(std::launch::async, [this](const std::string &path){return config_open(this->conf, path);}, conf_path);
}

void Aero::config_open_await() {
	Settings conf_settings = future_config_open.get();
	settings = conf_settings;
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
	// This is weird and might be slow
	static int type = static_cast<int>(aero.settings.rans.solver_type == "implicit");
	static const char* types[] = { "explicit", "implicit" };
	ImGui::RadioButton("Explicit", &type, 0); ImGui::SameLine();
	ImGui::RadioButton("Implicit", &type, 1);
	if (aero.settings.rans.solver_type != types[type]) {
		aero.settings.rans.solver_type = types[type];
	}
	static int order = 0;
	static const bool orders[] = { false, true };
	ImGui::RadioButton("First", &order, 0); ImGui::SameLine();
	ImGui::RadioButton("Second", &order, 1);
	if (aero.settings.rans.second_order != orders[order]) {
		aero.settings.rans.second_order = orders[order];
	}
	ImGui::InputDouble("Relaxation", &aero.settings.rans.relaxation, 0.1f, 1.0f, "%.2f");

	ImGui::End();
}

void ButtonLayer::OnUIRender() {
	{
		ImGui::Begin("Buttons");
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

		if (ImGui::Button("Open Config", ImVec2(-1.0f, 0.0f)) && !aero.signal_status_busy && aero.signal_status_ready) {
			// TODO: Make this a file dialog !!
			aero.config_open_async("../../../../examples/conf.ini");
		};

		if (!aero.signal_status_ready && is_future_done(aero.future_solve)) {
			aero.solve_await();
		}

		if (!aero.signal_status_ready && is_future_done(aero.future_config_open)) {
			aero.config_open_await();
		}

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

FlexGUI::Application* CreateApplication(int argc, char** argv, Aero& aero)
{
	FlexGUI::ApplicationSpecification spec;
	spec.Name = "GUI Example";
	spec.Width = 1600;
	spec.Height = 900;

	FlexGUI::Application* app = new FlexGUI::Application(spec);
	app->PushLayer(std::make_shared<ButtonLayer>(aero));
	app->PushLayer(std::make_shared<GraphLayer>(aero));
	app->PushLayer(std::make_shared<RansLayer>(aero));

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
		tiny::config conf;

		// Temporary
		rans::Settings data = rans::parse(argc, argv, __FILE__);
		if (data.read_failure == 1) {
			std::cout << "\033[0;31m";
			std::cout << "Error, unspecified error reading input file.";
			std::cout << "\033[0m\n" << std::endl;
			return 1;
		} else if (data.read_failure == 2) {
			std::cout << std::endl;
			return 1;
		}

		// Initialize modules with signal routing
		rans::Rans rans(gui);
		rans.settings = data;

		vlm::VLM vlm(gui);

		// Initialize main application with the modules
		Aero aero(conf, rans, vlm, gui);

		while (g_ApplicationRunning)
		{
			FlexGUI::Application* app = CreateApplication(argc, argv, aero);
			app->Run();
			delete app;
		}

		// rans.input();
		// rans.solve();

		return 0;
	}
}
