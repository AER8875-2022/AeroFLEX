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

template<typename T>
class VectorMutex {
	public:
		std::vector<T> m_vec;
		std::mutex m_mutex;

		VectorMutex() = default;
		VectorMutex(int size) : m_vec(size) {}

		T* data() {
			std::lock_guard<std::mutex> lock(m_mutex);
			if (m_vec.empty()) {
				return nullptr;
			}
			return &m_vec[0];
		}

		T& operator[](int i) {
			std::lock_guard<std::mutex> lock(m_mutex);
			return m_vec[i];
		}

		void push_back(const T& val) {
			std::lock_guard<std::mutex> lock(m_mutex);
			m_vec.push_back(val);
		}

		int size() {
			std::lock_guard<std::mutex> lock(m_mutex);
			return m_vec.size();
		}

		void resize(int size) {
			std::lock_guard<std::mutex> lock(m_mutex);
			m_vec.resize(size);
		}

		void clear() {
			std::lock_guard<std::mutex> lock(m_mutex);
			m_vec.clear();
		}
};

struct GUIHandler {
	SignalHandler signal;
	EventHandler event;
	spsc_queue<std::string> msg{8};
};
