#pragma once

#include <atomic>
#include <memory>
#include <string>

#include "spsc_queue.hpp"

struct SignalHandler {
	std::atomic<bool> stop = false;
	std::atomic<bool> pause = false;
};

struct EventHandler {
	// Rans
	// TODO: maybe move these into module specific structs
	std::atomic<bool> rans_preprocess = false;
	std::atomic<bool> rans_solve = false;
	std::atomic<bool> rans_postprocess = false;

	void reset() {
		rans_preprocess = false;
		rans_solve = false;
		rans_postprocess = false;
	}
};

struct GUIHandler {
	SignalHandler signal;
	EventHandler event;
	spsc_queue<std::string> msg{8};
};
