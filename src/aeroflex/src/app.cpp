#include "EntryPoint.h"

#include "imgui.h"
#include "implot.h"
#include "tinyconfig.hpp"
#include "parser.hpp"

#include "Application.h"
#include "FileDialog.hpp"

#include <string>
#include <cstring>
#include <future>
#include <thread>
#include <mutex>
#include <ctime>

#include "common_aeroflex.hpp"
#include "database/database.hpp"
#include "vlm/input.hpp"

#include <rans/rans.h>
#include <vlm/vlm.hpp>
#include <structure/structure.hpp>
#include <aero/aeroelasticity.hpp>


template <class T>
bool is_future_done(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

static void HelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
		ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

struct Settings {
	rans::Settings rans;
	vlm::Settings vlm;
    structure::Settings structure;
	aero::Settings aero;
};

enum class AppDialogAction {
	None,
	ConfigOpen,
	ConfigSave,
	DatabaseOpen,
	VlmMeshOpen,
};

class App {
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
        structure::Structure &structure;
		aero::Aero &aero;

		Settings settings;

		std::future<void> future_solve;
		std::future<std::optional<Settings>> future_config_open;
		std::future<bool> future_config_save;
		std::string outfile;

		// Signal routing
		bool signal_status_busy = false;
		bool signal_status_ready = true;
		GUIHandler &gui;

		// Dialogs
		AppDialogAction dialog_action = AppDialogAction::None;
		FlexGUI::FileDialog dialog;
		char path_buf[256];
		App(rans::Rans &rans, vlm::VLM &vlm, structure::Structure &structure, aero::Aero &aero, GUIHandler &gui);
};

struct GeoLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	GeoLayer(App &app) : app(app) {};
};

struct StructureLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	StructureLayer(App &app) : app(app) {};
};

struct RansLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	RansLayer(App &app) : app(app) {};
};

struct VlmLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	VlmLayer(App &app) : app(app) {};
};

struct ButtonLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	ButtonLayer(App &app) : app(app) {};
};

struct RansGraphLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	RansGraphLayer(App &app) : app(app) {};
};

struct StructureGraphLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	StructureGraphLayer(App &app) : app(app) {};
};

struct VlmGraphLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	VlmGraphLayer(App &app) : app(app) {};
};

struct AeroLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	AeroLayer(App &app) : app(app) {};
};

struct CpLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	CpLayer(App &app) : app(app) {};
};

struct ConsoleLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	ConsoleLayer(App &app) : app(app) {};
};

struct DialogLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	App &app;
	DialogLayer(App &app) : app(app) {};
};

/*
  #####  ####### #       #     # #######
 #     # #     # #       #     # #
 #       #     # #       #     # #
  #####  #     # #       #     # #####
       # #     # #        #   #  #
 #     # #     # #         # #   #
  #####  ####### #######    #    #######
*/
// =================================================================================================


void solve(rans::Rans &rans, vlm::VLM &vlm, structure::Structure &structure, aero::Aero &aero) {

	database::table table;



	if (!vlm.settings.sim.get_databaseFormat().compare("NONE")) {
		vlm.database.airfoils["naca0012q"];
		vlm.database.airfoils["naca0012q"].alpha = rans.alphas;

		for (auto& [airfoil, db] : vlm.database.airfoils) {
			rans.solve_airfoil(airfoil, db);
		}
	}
	else if (!vlm.settings.sim.get_databaseFormat().compare("FILE")) {
		vlm.database.importAirfoils(vlm.settings.io.databaseFile);
	}

	vlm.initialize();
	vlm.database.importLocations(vlm.settings.io.locationFile); // Temporary
	vlm.solve();

	structure.input();
	structure.solve();
}

// =================================================================================================

std::optional<Settings> config_open(const std::string &conf_path) {
	Settings settings;

	tiny::config io;
	bool success = io.read(conf_path);
	if (!success) return std::nullopt;

	settings.rans.import_config_file(io);
	settings.vlm.import_config_file(io);
	settings.structure.import_config_file(io);

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

	settings.rans.export_config_file(io);
	settings.vlm.export_config_file(io);
	settings.structure.export_config_file(io);

	return io.write(conf_path);
}

App::App(rans::Rans &rans, vlm::VLM &vlm, structure::Structure &structure, aero::Aero &aero, GUIHandler &gui) : rans(rans), vlm(vlm), structure(structure), aero(aero), gui(gui) {
	settings.rans.bcs["farfield"];
	settings.rans.bcs["farfield"].bc_type = "farfield";
	settings.rans.bcs["wall"];
	settings.rans.bcs["wall"].bc_type = "slip-wall";
}

void App::solve_async() {
	signal_status_ready = false;
	signal_status_busy = true;
	gui.msg.push("-- Starting simulation --");

	rans.settings = settings.rans;
	vlm.settings = settings.vlm;
	structure.settings = settings.structure;
	aero.settings = settings.aero;

	future_solve = std::async(std::launch::async,
	[&](){
		try {
			solve(rans, vlm, structure, aero);
		} catch (std::exception &e) {
			gui.msg.push(e.what());
		}
	});
}

void App::solve_await() {
	future_solve.get();
	if (gui.signal.stop) {
		gui.msg.push("-- Simulation stopped --");
	} else {
		gui.msg.push("-- Simulation done --");
	}
	signal_status_ready = true;
	signal_status_busy = false;
	gui.signal.stop = false;
}

void App::config_open_async(const std::string &conf_path) {
	signal_status_ready = false;
	signal_status_busy = true;
	outfile = conf_path;
	gui.msg.push("Starting parsing: " + conf_path);
	future_config_open = std::async(std::launch::async,
	[this](const std::string &path){
		std::optional<Settings> new_settings_op{};
		try {
			new_settings_op = config_open(path);
		} catch (std::exception &e) {
			gui.msg.push(e.what());
		}
		return new_settings_op;
	}, conf_path);
}

void App::config_open_await() {
	auto settings_op = future_config_open.get();
	if (settings_op.has_value()) {
		gui.msg.push("Config loaded.");
		settings = settings_op.value();
	} else {
		gui.msg.push("Config load failed.");
	}
	signal_status_ready = true;
	signal_status_busy = false;
}

void App::config_save_async(const std::string &conf_path) {
	signal_status_ready = false;
	signal_status_busy = true;
	gui.msg.push("Starting save: " + conf_path);
	future_config_save = std::async(std::launch::async,
	[this](const std::string &path){
		bool success = false;
		try {
			success = config_save(path, this->settings);
		} catch (std::exception &e) {
			gui.msg.push(e.what());
		}
		return success;
	}, conf_path);
}

void App::config_save_await() {
	bool success = future_config_save.get();
	if (success) {
		gui.msg.push("Config saved.");
	} else {
		gui.msg.push("Config save failed.");
	}
	signal_status_ready = true;
	signal_status_busy = false;
}

inline void Combo(std::vector<std::string> &vec, int &index, const char* label) {
	if (ImGui::BeginCombo(label, vec[index].c_str())) {
		for (int i = 0; i < vec.size(); i++) {
			bool is_selected = (index == i);
			if (ImGui::Selectable(vec[i].c_str(), is_selected)) index = i;
			if (is_selected) ImGui::SetItemDefaultFocus();
		}
		ImGui::EndCombo();
	}
}


void StructureLayer::OnUIRender() {
	ImGui::Begin("Structure");

	ImGui::Separator();
	ImGui::Text("Parameters");
	ImGui::InputDouble("Tolerance:", &app.settings.structure.Tolerance, 0.01f, 1.0f, "%.2e");

	ImGui::InputInt("N_step:", &app.settings.structure.N_step);

	ImGui::InputDouble("Damping", &app.settings.structure.Damping, 0.01f, 1.0f, "%.4f");

	ImGui::Separator();
	ImGui::Text("Options");
	if (ImGui::RadioButton("NONLINEAR", app.settings.structure.Solve_type == 0)) app.settings.structure.Solve_type = 0;
	ImGui::SameLine();
	if (ImGui::RadioButton("LINEAR", app.settings.structure.Solve_type == 1)) app.settings.structure.Solve_type = 1;
	ImGui::End();
}

void RansLayer::OnUIRender() {
	ImGui::Begin("RANS");

	ImGui::Separator();
	ImGui::Text("Gas");
	ImGui::InputDouble("gamma", &app.settings.rans.g.gamma, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Ratio of specific heats");
	ImGui::InputDouble("R", &app.settings.rans.g.R, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Specific gas constant");

	ImGui::Separator();
	ImGui::Text("Farfield");
	ImGui::InputDouble("Mach", &app.settings.rans.bcs["farfield"].vars_far.mach, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Mach number");
	ImGui::InputDouble("Temperature", &app.settings.rans.bcs["farfield"].vars_far.T, 0.01f, 1.0f, "%.4f");
	ImGui::InputDouble("Pressure", &app.settings.rans.bcs["farfield"].vars_far.p, 0.01f, 1.0f, "%.4f");

	ImGui::Separator();
	ImGui::Text("Alphas");
	ImGui::InputDouble("Alpha Start", &app.settings.rans.alpha_start, 0.1f, 1.0f, "%.1f");
	ImGui::SameLine(); HelpMarker("Alpha linspace start");
	ImGui::InputDouble("Alpha End", &app.settings.rans.alpha_end, 0.1f, 1.0f, "%.1f");
	ImGui::SameLine(); HelpMarker("Alpha linspace end");
	ImGui::InputDouble("Alpha Step", &app.settings.rans.alpha_step, 0.1f, 1.0f, "%.1f");
	ImGui::SameLine(); HelpMarker("Alpha linspace step");

	ImGui::Separator();
	ImGui::Text("Solver");

	ImGui::Checkbox("Second Order", &app.settings.rans.second_order);

	Combo(app.settings.rans.solver_options, app.settings.rans.solver, "Type");
	Combo(app.settings.rans.gradient_options, app.settings.rans.gradient, "Gradient");
	Combo(app.settings.rans.viscosity_options, app.settings.rans.viscosity, "Viscosity");

	ImGui::InputDouble("Relaxation", &app.settings.rans.relaxation, 0.1f, 1.0f, "%.2f");
	ImGui::InputDouble("Tolerance", &app.settings.rans.tolerance, 0.0f, 0.0f, "%e");
	ImGui::InputDouble("Start CFL", &app.settings.rans.start_cfl, 0.1f, 1.0f, "%.1f");
	ImGui::InputDouble("Slope CFL", &app.settings.rans.slope_cfl, 0.1f, 1.0f, "%.1f");
	ImGui::InputDouble("Limiter K", &app.settings.rans.limiter_k, 0.1f, 1.0f, "%.1f");
	ImGui::SameLine(); HelpMarker("Venkatakrishnan limiter parameter");
	ImGui::InputDouble("Max CFL", &app.settings.rans.max_cfl, 0.1f, 1.0f, "%.1f");
	ImGui::SliderInt("RHS Iterations", &app.settings.rans.rhs_iterations, 1, 10);
	ImGui::SameLine(); HelpMarker("Number of iterations for the RHS evaluation");
	ImGui::InputInt("Max Iterations", &app.settings.rans.max_iterations);

	ImGui::End();
}

void VlmLayer::OnUIRender() {
	ImGui::Begin("VLM");

	ImGui::Separator();
	ImGui::Text("Case Parameters");
	ImGui::InputDouble("AoA", &app.settings.vlm.sim.aoa, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Geometric angle of attack in degrees");
	ImGui::InputDouble("Sideslip", &app.settings.vlm.sim.sideslip, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Geometric angle of side slip in degrees");
	ImGui::InputDouble("V inf", &app.settings.vlm.sim.vinf, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Free stream magnitude velocity");
	ImGui::InputDouble("rho", &app.settings.vlm.sim.rho, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Density of the fluid");
	ImGui::InputDouble("cref", &app.settings.vlm.sim.cref, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Reference chord length");
	ImGui::InputDouble("sref", &app.settings.vlm.sim.sref, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Reference surface area");
	ImGui::InputDouble("coreRadius", &app.settings.vlm.sim.coreRadius, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Viscous relaxation value applied on the vortex filament kernel");
	ImGui::InputDouble("X ref", &app.settings.vlm.sim.x0, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("X component of origin to which the x and z moment are computed");
	ImGui::InputDouble("Y ref", &app.settings.vlm.sim.y0, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Y component of origin to which the x and z moment are computed");
	ImGui::InputDouble("Z ref", &app.settings.vlm.sim.z0, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Z component of origin to which the x and z moment are computed");

	Combo(app.settings.vlm.sim.databaseFormat_options, app.settings.vlm.sim.databaseFormat, "Db Format");
	if (app.settings.vlm.sim.get_databaseFormat() == "FILE") {
		ImGui::InputText("", app.settings.vlm.io.databaseFile.data(), app.settings.vlm.io.databaseFile.size(), ImGuiInputTextFlags_ReadOnly);
		ImGui::SameLine();
		if (ImGui::Button("...")) {
			app.dialog.file_dialog_open = true;
			app.dialog.type = FlexGUI::FileDialogType::OpenFile;
			app.dialog_action = AppDialogAction::DatabaseOpen;
		}
	}

	ImGui::Separator();
	ImGui::Text("Mesh");
	ImGui::InputText("", app.settings.vlm.io.meshFile.data(), app.settings.vlm.io.meshFile.size(), ImGuiInputTextFlags_ReadOnly);
	ImGui::SameLine();
	if (ImGui::Button("...")) {
		app.dialog.file_dialog_open = true;
		app.dialog.type = FlexGUI::FileDialogType::OpenFile;
		app.dialog_action = AppDialogAction::VlmMeshOpen;
	}

	ImGui::Separator();
	ImGui::Text("Solver");
	Combo(app.settings.vlm.solver.timeDomain_options, app.settings.vlm.solver.timeDomain, "Time Domain");
	Combo(app.settings.vlm.solver.type_options, app.settings.vlm.solver.type, "Type");
	Combo(app.settings.vlm.solver.linearSolver_options, app.settings.vlm.solver.linearSolver, "Linear solver");
	ImGui::InputDouble("Tolerance", &app.settings.vlm.solver.tolerance, 0.01f, 1.0f, "%e");
	ImGui::InputDouble("Relaxation", &app.settings.vlm.solver.relaxation, 0.01f, 1.0f, "%.4f");
	ImGui::InputInt("Max Iterations", &app.settings.vlm.solver.max_iter);

	ImGui::End();
};

void AeroLayer::OnUIRender() {
	ImGui::Begin("Aero");
	ImGui::Separator();
	ImGui::InputDouble("Tolerance", &app.settings.aero.tolerance, 0.01f, 1.0f, "%e");
	ImGui::End();
};

void DialogLayer::OnUIRender() {
	app.dialog.Show(app.path_buf);
	if (app.dialog.ready) {
		std::string path(app.path_buf);
		app.dialog.ready = false;
		if (app.dialog_action == AppDialogAction::ConfigOpen) {
			app.config_open_async(path);
		} else if (app.dialog_action == AppDialogAction::ConfigSave) {
			app.config_save_async(path);
		} else if (app.dialog_action == AppDialogAction::DatabaseOpen) {
			app.settings.vlm.io.databaseFile = path;
		} else if (app.dialog_action == AppDialogAction::VlmMeshOpen) {
			app.settings.vlm.io.meshFile = path;
		}
		strcpy(app.path_buf, "");
	}
}

void ButtonLayer::OnUIRender() {
	{
		ImGui::Begin("Buttons");
		ImGui::Text("Simulation");
		if (ImGui::Button("Start Simulation", ImVec2(-1.0f, 0.0f)) && !app.signal_status_busy && app.signal_status_ready) {
			app.solve_async();
		};

		if (ImGui::Button("Stop Simulation", ImVec2(-1.0f, 0.0f)) && app.signal_status_busy) {
			app.gui.signal.stop = true;
			app.gui.msg.push("Stopping simulation");
		};

		if (ImGui::Button("Pause", ImVec2(-1.0f, 0.0f)) && app.signal_status_busy) {
			app.gui.signal.pause = true;
			app.gui.msg.push("Pausing simulation");
		};

		if (ImGui::Button("Resume", ImVec2(-1.0f, 0.0f)) && app.signal_status_busy) {
			app.gui.signal.pause = false;
		};

		ImGui::Separator();
		ImGui::Text("Config");
		if (ImGui::Button("Open", ImVec2(-1.0f, 0.0f)) && !app.signal_status_busy && app.signal_status_ready) {
			app.dialog.file_dialog_open = true;
			app.dialog.type = FlexGUI::FileDialogType::OpenFile;
			app.dialog_action = AppDialogAction::ConfigOpen;
		};

		if (ImGui::Button("Save", ImVec2(-1.0f, 0.0f))) {
			app.dialog.file_dialog_open = true;
			app.dialog.type = FlexGUI::FileDialogType::SaveFile;
			app.dialog_action = AppDialogAction::ConfigSave;
		}

		if (!app.signal_status_ready && is_future_done(app.future_solve)) {
			app.solve_await();
		};

		if (!app.signal_status_ready && is_future_done(app.future_config_open)) {
			app.config_open_await();
		};

		if (!app.signal_status_ready && is_future_done(app.future_config_save)) {
			app.config_save_await();
		};
		ImGui::End();

		// ImGui::ShowDemoWindow();
		// ImPlot::ShowDemoWindow();
	};
};

void RansGraphLayer::OnUIRender() {
	{
		ImGui::Begin("Rans-Convergence");
		ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.0f, 0.0f));
		static ImPlotAxisFlags xflags = ImPlotAxisFlags_None;
		static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_RangeFit;

		if (ImPlot::BeginPlot("Convergence", ImVec2(-1,-1))) {
			ImPlot::SetupAxes("Iterations","Residual",xflags,yflags);
			ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, 0, app.rans.iters);
			ImPlot::SetupAxisZoomConstraints(ImAxis_X1, 11.0, app.rans.iters);
			ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

			if (!app.signal_status_ready) {
				ImPlot::SetupAxisLimits(ImAxis_X1, 0, app.rans.iters, ImPlotCond_Always);
			}
			ImPlot::PlotLine("L2 residual", app.rans.residuals.data(), app.rans.iters);
			ImPlot::EndPlot();
		}

		ImGui::End();
	}
};

void StructureGraphLayer::OnUIRender() {
	{
		ImGui::Begin("Structure-Convergence");
		ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.0f, 0.0f));
		static ImPlotAxisFlags xflags = ImPlotAxisFlags_None;
		static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_RangeFit;

		if (ImPlot::BeginPlot("Convergence", ImVec2(-1,-1))) {
			ImPlot::SetupAxes("Iterations","Residual",xflags,yflags);
			ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, 0, app.structure.iters);
			ImPlot::SetupAxisZoomConstraints(ImAxis_X1, 11.0, app.structure.iters);
			ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

			if (!app.signal_status_ready) {
				ImPlot::SetupAxisLimits(ImAxis_X1, 0, app.vlm.iters, ImPlotCond_Always);
			}
			ImPlot::PlotLine("L2 residual", app.structure.residuals.data(), app.structure.iters);
			ImPlot::EndPlot();
		}

		ImGui::End();
	}
}

void VlmGraphLayer::OnUIRender() {
	{
		ImGui::Begin("Vlm-Convergence");
		ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.0f, 0.0f));
		static ImPlotAxisFlags xflags = ImPlotAxisFlags_None;
		static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_RangeFit;

		if (ImPlot::BeginPlot("Convergence", ImVec2(-1,-1))) {
			ImPlot::SetupAxes("Iterations","Residual",xflags,yflags);
			ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, 0, app.vlm.iters);
			ImPlot::SetupAxisZoomConstraints(ImAxis_X1, 11.0, app.vlm.iters);
			ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

			if (!app.signal_status_ready) {
				ImPlot::SetupAxisLimits(ImAxis_X1, 0, app.vlm.iters, ImPlotCond_Always);
			}
			ImPlot::PlotLine("L2 residual", app.vlm.residuals.data(), app.vlm.iters);
			ImPlot::EndPlot();
		}

		ImGui::End();
	}
};

void CpLayer::OnUIRender() {
	ImGui::Begin("Rans-Cp");
	auto size = ImGui::GetWindowSize();
	ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.2f, 0.2f));
	static ImPlotAxisFlags xflag = ImPlotAxisFlags_AutoFit;
	static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_Invert;

	if (ImPlot::BeginPlot("Cp Profile", ImVec2(-1, (int)(size.y / 2.2f)))) {
		std::scoped_lock lock(app.rans.profile.m_mutex);
		ImPlot::SetupAxes("x","Cp",xflag,yflags);
		ImPlot::PlotLine("Cp", app.rans.profile.x.data(), app.rans.profile.cp.data(), app.rans.profile.x.size());
		ImPlot::EndPlot();
	}

	if (ImPlot::BeginPlot("Cp Airfoil", ImVec2(-1, (int)(size.y / 2.2f)))) {
		const double pad = 1.1;
		std::scoped_lock lock(app.rans.profile.m_mutex);
		ImPlot::SetupAxes("x","", ImPlotAxisFlags_None,ImPlotAxisFlags_None);
		ImPlot::SetupAxisLimits(ImAxis_X1, app.rans.profile.xmin * pad, app.rans.profile.xmax * pad, ImPlotCond_Always);
		ImPlot::SetupAxisLimits(ImAxis_Y1, app.rans.profile.ymin * pad, app.rans.profile.ymax * pad, ImPlotCond_Always);
		ImPlot::PlotLine("Airfoil", app.rans.profile.x.data(), app.rans.profile.y.data(), app.rans.profile.x.size());

		if (app.rans.profile.filled) {
			ImPlot::PushPlotClipRect();
			for (int i = 0; i < app.rans.profile.x.size(); i++) {
				ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint(app.rans.profile.x[i], app.rans.profile.y[i]));
				ImVec2 p2 = ImPlot::PlotToPixels(ImPlotPoint(app.rans.profile.cp_airfoil_x[i], app.rans.profile.cp_airfoil_y[i]));
				if (app.rans.profile.cp_positive[i] == 1) {
					ImPlot::GetPlotDrawList()->AddLine(p1, p2, IM_COL32(0,255,0,255));
				} else {
					ImPlot::GetPlotDrawList()->AddLine(p1, p2, IM_COL32(255,0,0,255));
				}
			}
			ImPlot::PopPlotClipRect();
		}
		ImPlot::EndPlot();
	}

	ImGui::End();
};

class ConsoleLog
{
	public:
		void log(const char* msg) {
			time_t ttime = time(0);
			tm *t = localtime(&ttime);
			add_entry("[%02d:%02d:%02d] : %s\n", t->tm_hour, t->tm_min, t->tm_sec, msg);
		}

		void draw(const char* title, bool* p_open = NULL) {
			if (!ImGui::Begin(title, p_open)) {
				ImGui::End();
				return;
			}

			// Main window
			bool clear = ImGui::Button("Clear");
			ImGui::SameLine();
			bool copy = ImGui::Button("Copy");
			ImGui::SameLine();
			Filter.Draw("Filter", -100.0f);

			ImGui::Separator();

			if (ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar))
			{
				if (clear)
					clear_console();
				if (copy)
					ImGui::LogToClipboard();

				ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
				const char* buf = Buf.begin();
				const char* buf_end = Buf.end();
				if (Filter.IsActive()) {
					for (int line_no = 0; line_no < LineOffsets.Size; line_no++)
					{
						const char* line_start = buf + LineOffsets[line_no];
						const char* line_end = (line_no + 1 < LineOffsets.Size) ? (buf + LineOffsets[line_no + 1] - 1) : buf_end;
						if (Filter.PassFilter(line_start, line_end))
							ImGui::TextUnformatted(line_start, line_end);
					}
				} else {
					ImGuiListClipper clipper;
					clipper.Begin(LineOffsets.Size);
					while (clipper.Step()) {
						for (int line_no = clipper.DisplayStart; line_no < clipper.DisplayEnd; line_no++) {
							const char* line_start = buf + LineOffsets[line_no];
							const char* line_end = (line_no + 1 < LineOffsets.Size) ? (buf + LineOffsets[line_no + 1] - 1) : buf_end;
							ImGui::TextUnformatted(line_start, line_end);
						}
					}
					clipper.End();
				}
				ImGui::PopStyleVar();

				if (AutoScroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
					ImGui::SetScrollHereY(1.0f);
			}
			ImGui::EndChild();
			ImGui::End();
		}

		void add_entry(const char* fmt, ...) IM_FMTARGS(2)
		{
			int old_size = Buf.size();
			va_list args;
			va_start(args, fmt);
			Buf.appendfv(fmt, args);
			va_end(args);
			for (int new_size = Buf.size(); old_size < new_size; old_size++)
				if (Buf[old_size] == '\n')
					LineOffsets.push_back(old_size + 1);
		};

		ConsoleLog() {
			log(""); // Hack
			log("AeroFLEX v0.1");
		}

	private:
		ImGuiTextBuffer Buf;
		ImGuiTextFilter Filter;
		ImVector<int> LineOffsets; // Index to lines offset. We maintain this with AddLog() calls.
    	const bool AutoScroll = true;  // Keep scrolling if already at the bottom.

		void clear_console() {
			Buf.clear();
			LineOffsets.clear();
			LineOffsets.push_back(0);
		};
};

void ConsoleLayer::OnUIRender() {
    static ConsoleLog log;

    ImGui::Begin("Console", NULL);
    std::optional<std::string> txt = app.gui.msg.pop();

    if (txt.has_value()) {
        log.log(txt.value().c_str());
    }
    ImGui::End();
    log.draw("Console", NULL);
}

FlexGUI::Application* CreateApplication(int argc, char** argv, App& app)
{
	FlexGUI::ApplicationSpecification spec;
	spec.Name = "AEROFLEX";
	spec.Width = 1600;
	spec.Height = 900;

	FlexGUI::Application* application = new FlexGUI::Application(spec);
	application->SetTheme(FlexGUI::Theme::Light);
	application->PushLayer(std::make_shared<ButtonLayer>(app));
	application->PushLayer(std::make_shared<RansGraphLayer>(app));
	application->PushLayer(std::make_shared<VlmGraphLayer>(app));
	application->PushLayer(std::make_shared<StructureGraphLayer>(app));
	application->PushLayer(std::make_shared<CpLayer>(app));
	application->PushLayer(std::make_shared<StructureLayer>(app));
	application->PushLayer(std::make_shared<RansLayer>(app));
	application->PushLayer(std::make_shared<VlmLayer>(app));
	application->PushLayer(std::make_shared<AeroLayer>(app));
	application->PushLayer(std::make_shared<ConsoleLayer>(app));
	application->PushLayer(std::make_shared<DialogLayer>(app));

	application->SetMenubarCallback([application]()
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Exit"))
			{
				application->Close();
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Options"))
		{
			if (ImGui::BeginMenu("Theme")) {
				if (ImGui::MenuItem("Dark")) {
					application->SetTheme(FlexGUI::Theme::Dark);
				}
				if (ImGui::MenuItem("Light")) {
					application->SetTheme(FlexGUI::Theme::Light);
				}
				ImGui::EndMenu();
			}
			ImGui::EndMenu();
		}
		ImGui::SameLine(ImGui::GetWindowWidth() - 70);
		ImGui::Text("FPS: %.0f", 1.0f / ImGui::GetIO().DeltaTime);
	});
	return application;
}

void cmdl_parser_configure(cmd_line_parser::parser& parser) {
    parser.add("config",                 // name
               "Configuration file",  // description
               "-i",                   // shorthand
               false,                   // required argument
               false                   // is boolean option
    );
}

namespace FlexGUI {
	int Main(int argc, char** argv) {
		cmd_line_parser::parser parser(argc, argv);
		cmdl_parser_configure(parser);

		bool success = parser.parse();
		if (!success) return 1;
		std::string config_file = parser.get<std::string>("config");

		GUIHandler gui;

		// Initialize modules with signal routing
		rans::Rans rans(gui);
		vlm::VLM vlm(gui);
		structure::Structure structure(gui);
		aero::Aero aero(gui, vlm, structure);

		// Initialize main application with the modules

		App app(rans, vlm, structure, aero, gui);

		if (config_file == "") {
			while (g_ApplicationRunning) {
				FlexGUI::Application* application = CreateApplication(argc, argv, app);
				application->Run();
				delete application;
			}
		} else {
			std::cout << "CLI mode" << std::endl;
			std::optional<Settings> settings_opt;
			try {
				settings_opt = config_open(parser.get<std::string>("config"));
			} catch (std::exception &e) {
				std::cout << e.what() << std::endl;
				return 1;
			}
			if (!settings_opt.has_value()) return 1;
			app.settings = settings_opt.value();
			app.solve_async();
			while (g_ApplicationRunning) {
				if (is_future_done(app.future_solve)) {
					app.solve_await();
					g_ApplicationRunning = false;
				};
				auto txt = app.gui.msg.pop();
				while (txt.has_value()) {
					std::cout << txt.value() << std::endl;
					txt = app.gui.msg.pop();
				}

				std::this_thread::sleep_for(std::chrono::milliseconds(200));
			}
		}

		return 0;
	}
}
