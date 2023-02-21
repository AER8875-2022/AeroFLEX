#pragma once

#include <atomic>
#include <memory>

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
class LockFreeQueue
{
private:
    struct node
    {
        std::shared_ptr<T> data;
        std::atomic<node*> next;
        node(T const& data_):
            data(std::make_shared<T>(data_))
        {}
    };
    std::atomic<node*> head;
    std::atomic<node*> tail;
public:
    void push(T const& data)
    {
        std::atomic<node*> const new_node=new node(data);
        node* old_tail = tail.load();
        while(!old_tail->next.compare_exchange_weak(nullptr, new_node)){
          node* old_tail = tail.load();
        }
        tail.compare_exchange_weak(old_tail, new_node);
    }
    std::shared_ptr<T> pop()
    {
        node* old_head=head.load();
        while(old_head &&
            !head.compare_exchange_weak(old_head,old_head->next)){
            old_head=head.load();
        }
        return old_head ? old_head->data : std::shared_ptr<T>();
    }
};

struct GUIHandler {
	SignalHandler signal;
	EventHandler event;
	LockFreeQueue<std::string> message;
};
