#pragma once

#include <atomic>

struct SignalHandler {
	std::atomic<bool> stop = false;
	std::atomic<bool> pause = false;
};