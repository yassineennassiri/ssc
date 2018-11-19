#ifndef _SINGLETON_HOLDER_H_
#define _SINGLETON_HOLDER_H_

#include <cstdlib>
#include <memory>
#include <stdexcept>

template <class T>
struct CreateUsingUniquePtr
{
	static T * Create() {
		std::unique_ptr<T> tmp(new T);
		return tmp->get();
	}
};

template <class T>
struct CreateUsingNew
{
	static T * Create() {
		return new T;
	}
	static void Destroy(T* p) {
		delete p;
	}
};


/**
For now implement SingletonHolder with only the CreationPolicy
Next steps would be to add LifetimePolicy and ThreadingModel
*/
template
<
	class T,
	template <class> class CreationPolicy = CreateUsingNew 
>
class SingletonHolder
{
public:
	typedef T ObjectType;

	static T& Instance() {
		if (!pInstance_) {
			MakeInstance();
		}
		return pInstance_;
	}

private:

	// Helpers
	static void MakeInstance() {
		if (!pInstance_) {
			if (destroyed_) {
				destroyed_ = false;
			}
			pInstance_ = CreationPolicy<T>::Create();
		}
	}
	static void DestroySingleton() {
		if (!destroyed_) {
			CreationPolicy<T>::Destroy(pInstance_);
			pInstance_ = nullptr;
			destroyed_ = true;
		}

	}

	// Protection
	SingletonHolder() {};
	~SingletonHolder() {};

	// Data
	static T * pInstance_;
	static bool destroyed_;
};

#endif