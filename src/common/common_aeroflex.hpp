#pragma once

#include <atomic>
#include <memory>
#include <string>
#include <mutex>
#include <vector>

#include "spsc_queue.hpp"

struct SignalHandler {
	std::atomic<bool> stop = false;
	std::atomic<bool> pause = false;
};

struct GUIHandler {
	SignalHandler signal;
	spsc_queue<std::string> msg{8};
};
