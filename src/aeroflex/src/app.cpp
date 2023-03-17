#include "EntryPoint.h"

#include "imgui.h"
#include "implot.h"
#include "tinyconfig.hpp"

#include "Application.h"
#include "FileDialog.hpp"

#include <string>
#include <cstring>
#include <future>
#include <thread>
#include <ctime>

#include "common_aeroflex.hpp"
#include "database/database.hpp"
#include "vlm/input.hpp"

#include <rans/rans.h>
#include <vlm/vlm.hpp>

template <class T>
bool is_future_done(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

const std::string bool_to_string(const bool b) {
    return b ? "true" : "false";
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
		GUIHandler &gui;
		Aero(rans::Rans &rans, vlm::VLM &vlm, GUIHandler &gui);
};

struct RansLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	RansLayer(Aero &aero) : aero(aero) {};
};

struct VlmLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	VlmLayer(Aero &aero) : aero(aero) {};
};

struct ButtonLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	ButtonLayer(Aero &aero) : aero(aero) {};
};

struct RansGraphLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	RansGraphLayer(Aero &aero) : aero(aero) {};
};

struct VlmGraphLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	VlmGraphLayer(Aero &aero) : aero(aero) {};
};

struct CpLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	CpLayer(Aero &aero) : aero(aero) {};
};

struct ConsoleLayer : public FlexGUI::Layer {
	virtual void OnUIRender() override;
	Aero &aero;
	ConsoleLayer(Aero &aero) : aero(aero) {};
};

// SOLVE =================================================================================================

void solve(rans::Rans &rans, vlm::VLM &vlm) {

	database::table table;
	table.airfoils["naca0012_coarse"];
	table.airfoils["naca0012_coarse"].alpha = {1.0, 5.0, 10.0};

	for (auto& [airfoil, db] : table.airfoils) {
		rans.solve_airfoil(airfoil, db);
	}

	vlm.input();
	vlm.database = table;
	vlm.database.importLocations(vlm.settings.io.locationFile);
	vlm.solve();

	// rans.input();
	// rans.solve();
}

// =================================================================================================

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

	settings.rans.set_solver_type(io.get<std::string>("rans-solver", "solver"));
	settings.rans.set_gradient_scheme(io.get<std::string>("rans-solver", "gradient"));
	settings.rans.set_viscosity_model(io.get<std::string>("rans-solver", "viscosity"));
	settings.rans.second_order = io.get<bool>("rans-solver", "second_order");
	settings.rans.relaxation = io.get<double>("rans-solver", "relaxation");
	settings.rans.start_cfl = io.get<double>("rans-solver", "start_cfl");
	settings.rans.slope_cfl = io.get<double>("rans-solver", "slope_cfl");
	settings.rans.max_cfl = io.get<double>("rans-solver", "max_cfl");
	settings.rans.tolerance = io.get<double>("rans-solver", "tolerance");
	settings.rans.rhs_iterations = io.get<int>("rans-solver", "rhs_iterations");
	settings.rans.max_iterations = io.get<int>("rans-solver", "max_iterations");
	settings.rans.limiter_k = io.get<double>("rans-solver", "limiter_k");
 
	for (int i = 0; i < io.how_many("rans-mesh"); i++) {
		settings.rans.meshes.push_back(io.get_i<std::string>("rans-mesh", "file", i));
	}

	// Parsing for vlm
	settings.vlm.import_config_file(io);

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
		std::cout << name << " | " << bc.bc_type << std::endl;
		io.config_vec["rans-bc"].push_back({
			{"type", bc.bc_type},
			{"name", name},
			{"T", std::to_string(bc.vars_far.T)},
			{"mach", std::to_string(bc.vars_far.mach)},
			{"angle", std::to_string(bc.vars_far.angle)},
			{"p", std::to_string(bc.vars_far.p)},
		});
	};

	io.config["rans-solver"]["solver"] = settings.rans.solver_type();
	io.config["rans-solver"]["gradient"] = settings.rans.gradient_scheme();
	io.config["rans-solver"]["viscosity"] = settings.rans.viscosity_model();

	io.config["rans-solver"]["second_order"] = bool_to_string(settings.rans.second_order);
	io.config["rans-solver"]["relaxation"] = bool_to_string(settings.rans.relaxation);
	io.config["rans-solver"]["start_cfl"] = std::to_string(settings.rans.start_cfl);
	io.config["rans-solver"]["slope_cfl"] = std::to_string(settings.rans.slope_cfl);
	io.config["rans-solver"]["max_cfl"] = std::to_string(settings.rans.max_cfl);
	io.config["rans-solver"]["tolerance"] = std::to_string(settings.rans.tolerance);
	io.config["rans-solver"]["rhs_iterations"] = std::to_string(settings.rans.rhs_iterations);
	io.config["rans-solver"]["max_iterations"] = std::to_string(settings.rans.max_iterations);
	io.config["rans-solver"]["limiter_k"] = std::to_string(settings.rans.limiter_k);

	io.config_vec["rans-mesh"] = {};
	for (auto &mesh : settings.rans.meshes) {
		io.config_vec["rans-mesh"].push_back({
			{"file", mesh},
		});
	};

	// Exporting for VLM
	settings.vlm.export_config_file(io);

	return io.write(conf_path);
}

Aero::Aero(rans::Rans &rans, vlm::VLM &vlm, GUIHandler &gui) : rans(rans), vlm(vlm), gui(gui) {
	settings.rans.bcs["farfield"];
	settings.rans.bcs["farfield"].bc_type = "farfield";
	settings.rans.bcs["wall"];
	settings.rans.bcs["wall"].bc_type = "slip-wall";

	// TEMPORARY !!!!
	settings.rans.meshes.push_back("../../../../examples/rans/airfoil_API2_coarse.msh");
	settings.rans.meshes.push_back("../../../../examples/rans/airfoil_API2_mid.msh");
	settings.rans.meshes.push_back("../../../../examples/rans/airfoil_API2_fine.msh");
}

void Aero::solve_async() {
	signal_status_ready = false;
	signal_status_busy = true;
	gui.msg.push("-- Starting simulation --");
	rans.settings = settings.rans;
	vlm.settings = settings.vlm;
	future_solve = std::async(std::launch::async, 
	[&](){
		try {
			solve(rans, vlm);
		} catch (std::exception &e) {
			gui.msg.push(e.what());
		}
	});
}

void Aero::solve_await() {
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

void Aero::config_open_async(const std::string &conf_path) {
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

void Aero::config_open_await() {
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

void Aero::config_save_async(const std::string &conf_path) {
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

void Aero::config_save_await() {
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

	ImGui::Checkbox("Second Order", &aero.settings.rans.second_order);

	Combo(aero.settings.rans.solver_options, aero.settings.rans.solver, "Type");
	Combo(aero.settings.rans.gradient_options, aero.settings.rans.gradient, "Gradient");
	Combo(aero.settings.rans.viscosity_options, aero.settings.rans.viscosity, "Viscosity");

	ImGui::InputDouble("Relaxation", &aero.settings.rans.relaxation, 0.1f, 1.0f, "%.2f");
	ImGui::InputDouble("Tolerance", &aero.settings.rans.tolerance, 0.0f, 0.0f, "%e");
	ImGui::InputDouble("Start CFL", &aero.settings.rans.start_cfl, 0.1f, 1.0f, "%.1f");
	ImGui::InputDouble("Slope CFL", &aero.settings.rans.slope_cfl, 0.1f, 1.0f, "%.1f");
	ImGui::InputDouble("Limiter K", &aero.settings.rans.limiter_k, 0.1f, 1.0f, "%.1f");
	ImGui::InputDouble("Max CFL", &aero.settings.rans.max_cfl, 0.1f, 1.0f, "%.1f");
	ImGui::SliderInt("RHS Iterations", &aero.settings.rans.rhs_iterations, 1, 10);
	ImGui::InputInt("Max Iterations", &aero.settings.rans.max_iterations);

	ImGui::End();
}

void VlmLayer::OnUIRender() {
	ImGui::Begin("VLM");

	ImGui::Separator();
	ImGui::Text("Case Parameters");
	ImGui::InputDouble("Sideslip", &aero.settings.vlm.sim.sideslip, 0.01f, 1.0f, "%.4f"); 
	ImGui::SameLine(); HelpMarker("Geometric angle of side slip in degrees");
	ImGui::InputDouble("V inf", &aero.settings.vlm.sim.vinf, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Free stream magnitude velocity");
	ImGui::InputDouble("rho", &aero.settings.vlm.sim.rho, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Density of the fluid");
	ImGui::InputDouble("cref", &aero.settings.vlm.sim.cref, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Reference chord length");
	ImGui::InputDouble("sref", &aero.settings.vlm.sim.sref, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Reference surface area");
	ImGui::InputDouble("coreRadius", &aero.settings.vlm.sim.coreRadius, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Viscous relaxation value applied on the vortex filament kernel");
	ImGui::InputDouble("X ref", &aero.settings.vlm.sim.x0, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("X component of origin to which the x and z moment are computed");
	ImGui::InputDouble("Y ref", &aero.settings.vlm.sim.y0, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Y component of origin to which the x and z moment are computed");
	ImGui::InputDouble("Z ref", &aero.settings.vlm.sim.z0, 0.01f, 1.0f, "%.4f");
	ImGui::SameLine(); HelpMarker("Z component of origin to which the x and z moment are computed");
	Combo(aero.settings.vlm.sim.databaseFormat_options, aero.settings.vlm.sim.databaseFormat, "Db Format");
	
	ImGui::Separator();
	ImGui::Text("Solver");
	Combo(aero.settings.vlm.solver.timeDomain_options, aero.settings.vlm.solver.timeDomain, "Time Domain");
	Combo(aero.settings.vlm.solver.type_options, aero.settings.vlm.solver.type, "Type");
	Combo(aero.settings.vlm.solver.linearSolver_options, aero.settings.vlm.solver.linearSolver, "Linear solver");
	ImGui::InputDouble("Tolerance", &aero.settings.vlm.solver.tolerance, 0.01f, 1.0f, "%e");
	ImGui::InputDouble("Relaxation", &aero.settings.vlm.solver.relaxation, 0.01f, 1.0f, "%.4f");
	ImGui::InputInt("Max Iterations", &aero.settings.vlm.solver.max_iter);

	ImGui::End();
}

void ButtonLayer::OnUIRender() {
	{
		static FlexGUI::FileDialog fd;
		ImGui::Begin("Buttons");
		ImGui::Text("Simulation");
		if (ImGui::Button("Start Simulation", ImVec2(-1.0f, 0.0f)) && !aero.signal_status_busy && aero.signal_status_ready) {
			aero.solve_async();
		};

		if (ImGui::Button("Stop Simulation", ImVec2(-1.0f, 0.0f)) && aero.signal_status_busy) {
			aero.gui.signal.stop = true;
			aero.gui.msg.push("Stopping simulation");
		};

		if (ImGui::Button("Pause", ImVec2(-1.0f, 0.0f)) && aero.signal_status_busy) {
			aero.gui.signal.pause = true;
			aero.gui.msg.push("Pausing simulation");
		};

		if (ImGui::Button("Resume", ImVec2(-1.0f, 0.0f)) && aero.signal_status_busy) {
			aero.gui.signal.pause = false;
		};

		ImGui::Separator();
		ImGui::Text("Config");
		if (ImGui::Button("Open", ImVec2(-1.0f, 0.0f)) && !aero.signal_status_busy && aero.signal_status_ready) {
			fd.file_dialog_open = true;
			fd.type = FlexGUI::FileDialogType::OpenFile;
		};

		if (ImGui::Button("Save", ImVec2(-1.0f, 0.0f))) {
			fd.file_dialog_open = true;
			fd.type = FlexGUI::FileDialogType::SaveFile;
		}

		if (!aero.signal_status_ready && is_future_done(aero.future_solve)) {
			aero.solve_await();
		};

		if (!aero.signal_status_ready && is_future_done(aero.future_config_open)) {
			aero.config_open_await();
		};

		if (!aero.signal_status_ready && is_future_done(aero.future_config_save)) {
			aero.config_save_await();
		};

		static char path_buf[256];

		fd.Show(path_buf);

		if (fd.ready) {
			std::string path(path_buf);
			fd.ready = false;
			// TODO: replace these with custom enums (eventually)
			if (fd.type == FlexGUI::FileDialogType::OpenFile) {
				aero.config_open_async(path);
			} else if (fd.type == FlexGUI::FileDialogType::SaveFile) {
				aero.config_save_async(path);
			}
			strcpy(path_buf, "");
		}
		ImGui::End();

		// ImGui::ShowDemoWindow();
		// ImPlot::ShowDemoWindow();
	};
};

void RansGraphLayer::OnUIRender() {
	{
		ImGui::Begin("Rans-Convergence");
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

void VlmGraphLayer::OnUIRender() {
	{
		ImGui::Begin("Vlm-Convergence");
		static ImPlotAxisFlags xflags = ImPlotAxisFlags_None;
		static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_RangeFit;
		const double xticks = 1;

		if (ImPlot::BeginPlot("Convergence", ImVec2(-1,400))) {
			ImPlot::SetupAxes("Iterations","Residual",xflags,yflags);
			ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, 0, aero.vlm.iters);
			ImPlot::SetupAxisZoomConstraints(ImAxis_X1, 11.0, aero.vlm.iters);
			ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

			if (!aero.signal_status_ready) {
				ImPlot::SetupAxisLimits(ImAxis_X1, 0, aero.vlm.iters, ImPlotCond_Always);
			}
			ImPlot::PlotLine("L2 residual", aero.vlm.residuals.data(), aero.vlm.iters);
			ImPlot::EndPlot();
		}

		ImGui::End();
	}
};

void CpLayer::OnUIRender() {
	{
		ImGui::Begin("Rans-Cp");
		static ImPlotAxisFlags xflags = ImPlotAxisFlags_None;
		static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit|ImPlotAxisFlags_RangeFit|ImPlotAxisFlags_Invert;
		const double xticks = 1;

		if (ImPlot::BeginPlot("Convergence", ImVec2(-1,400))) {
			ImPlot::SetupAxisLimits(ImAxis_X1, -0.1, 1.1, ImPlotCond_Always);
			ImPlot::SetupAxes("x","Cp",xflags,yflags);
			ImPlot::PlotLine("Cp", aero.rans.profile.x.data(), aero.rans.profile.cp.data(), aero.rans.profile.x.size());
			ImPlot::EndPlot();
		}

		ImGui::End();
	}
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
    std::optional<std::string> txt = aero.gui.msg.pop();

    if (txt.has_value()) {
        log.log(txt.value().c_str());
    }
    ImGui::End();
    log.draw("Console", NULL);
}

FlexGUI::Application* CreateApplication(int argc, char** argv, Aero& aero)
{
	FlexGUI::ApplicationSpecification spec;
	spec.Name = "AEROFLEX";
	spec.Width = 1600;
	spec.Height = 900;

	FlexGUI::Application* app = new FlexGUI::Application(spec);
	app->PushLayer(std::make_shared<ButtonLayer>(aero));
	app->PushLayer(std::make_shared<RansGraphLayer>(aero));
	app->PushLayer(std::make_shared<VlmGraphLayer>(aero));
	app->PushLayer(std::make_shared<CpLayer>(aero));
	app->PushLayer(std::make_shared<RansLayer>(aero));
	app->PushLayer(std::make_shared<VlmLayer>(aero));
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
	int Main(int argc, char** argv) {
		GUIHandler gui;

		// Initialize modules with signal routing
		rans::Rans rans(gui);
		vlm::VLM vlm(gui);

		// Initialize main application with the modules
		Aero aero(rans, vlm, gui);

		while (g_ApplicationRunning) {
			FlexGUI::Application* app = CreateApplication(argc, argv, aero);
			app->Run();
			delete app;
		}
		return 0;
	}
}
