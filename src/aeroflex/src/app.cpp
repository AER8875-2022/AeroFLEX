#include "Application.h"
#include "EntryPoint.h"
#include "imgui.h"
#include "implot.h"

#include "common_aeroflex.hpp"

#include <rans/rans.h>

// Temporary
#include <rans/parser.h>

#include <future>
#include <thread>

using namespace rans;

class Aero {
	public:
		void solve();
		void solve_async();

		// Modules
		Rans &rans;

		double x = 0.0;

		std::future<void> future_solve;

		// Signal routing
		bool signal_status_busy = false;
		bool signal_status_ready = true;
		GUIHandler &gui;
		Aero(Rans &rans, GUIHandler &gui) : rans(rans), gui(gui) {};
};

template <class T>
bool is_future_done(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

void Aero::solve_async() {
	signal_status_ready = false;
	signal_status_busy = true;
	future_solve = std::async(std::launch::async, [this](){solve();});
}

void Aero::solve() {
	rans.input();
	rans.solve();
}

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

		if (!aero.signal_status_ready && is_future_done(aero.future_solve)) {
			aero.future_solve.get();
			aero.signal_status_ready = true;
			aero.signal_status_busy = false;
			aero.gui.signal.stop = false;
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

		// Temporary
		Settings data = parse(argc, argv, __FILE__);
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
		Rans rans(gui);
		rans.data = data;

		// Initialize main application with the modules
		Aero aero(rans, gui);

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
