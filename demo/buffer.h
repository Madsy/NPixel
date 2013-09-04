#ifndef BUFFER_H_GUARD
#define BUFFER_H_GUARD
#include <SDL/SDL.h>
#include "myassert.h"
#ifdef DEBUG
#include <stdexcept>
#endif

template<typename T> struct Buffer2D {
	Buffer2D(unsigned int width,
	         unsigned int height,
	         T* buf = 0,
	         bool hw = false) :
		isHWBuffer(hw),
		w(width),
		h(height) {
		if(isHWBuffer) {
			data = 0;
		} else {
			data = new T[width * height];
			if(buf)
				memcpy(data, buf, sizeof(T) * width * height);
		}
	}

	Buffer2D(const Buffer2D<T>& b) {
		/*
		if(!isHWBuffer){
		  delete [] data;
		  data = 0;
		}
		*/
		isHWBuffer = b.isHWBuffer;
		w = b.w;
		h = b.h;

		if(isHWBuffer)
			data = 0;
		else {
			data = new T[w * h];
			if(b.data != 0)
				memcpy(data, b.data, w*h*sizeof(T));
		}
	}
	Buffer2D<T>& operator=(const Buffer2D<T>& b) {
		if(!isHWBuffer) {
			delete [] data;
			data = 0;
		}
		isHWBuffer = b.isHWBuffer;
		w = b.w;
		h = b.h;
		if(isHWBuffer)
			data = 0;
		else {
			data = new T[w * h];
			if(b.data != 0)
				memcpy(data, b.data, w*h*sizeof(T));
		}
		return *this;
	}
	Buffer2D() : isHWBuffer(false), w(0), h(0), data(0) {
		printf("Creating default buffer\n");
	}
	~Buffer2D() {
		if(!isHWBuffer) {
			delete[] data;
		}
		data = 0;
	}
	inline T& operator[](size_t index) {
		return data[index];
	}

	inline const T& operator[](size_t index) const {
		return data[index];
	}
	T* Lock() {
		if(isHWBuffer) {
			SDL_Surface* s = SDL_GetVideoSurface();
			SDL_LockSurface(s);
			return (T*)s->pixels;
		}
		return data;
	}
	void Unlock() {
		if(isHWBuffer) {
			SDL_Surface* s = SDL_GetVideoSurface();
			SDL_UnlockSurface(s);
		}
	}
	T* Ptr() {
		return data;
	}
	bool isHWBuffer;
	unsigned int w;
	unsigned int h;
	T* data;
};

template<typename T>
class Buffer1D {
public:
	Buffer1D(unsigned int capacity) : elemt_size(0), alloc_size(capacity), data(new T[capacity]) {}
	Buffer1D() : elemt_size(0), alloc_size(0), data(0) {}
	~Buffer1D() {
		delete[] data;
	}
	inline T& operator[](size_t index) {
#ifdef DEBUG
		if(index >= elemt_size) {
			std::runtime_error e("Buffer1D index out of bounds.");
			throw e;
		}
#endif
		return data[index];
	}
	inline const T& operator[](size_t index) const {
#ifdef DEBUG
		if(index >= elemt_size) {
			std::runtime_error e("Buffer1D index out of bounds.");
			throw e;
		}
#endif
		return data[index];
	}

	void Grow() {
		const unsigned int bucketSize = 1000;
		unsigned int new_alloc_size = alloc_size + bucketSize;
		T* newData = new T[new_alloc_size];
		memcpy(newData, data, sizeof(T) * elemt_size);
		alloc_size = new_alloc_size;
		delete [] data;
		data = newData;
	}
	void Grow(unsigned int extracap) {
		if(!extracap) return;
		unsigned int new_alloc_size = alloc_size + extracap;
		T* newData = new T[new_alloc_size];
		memcpy(newData, data, sizeof(T) * elemt_size);
		alloc_size = new_alloc_size;
		delete [] data;
		data = newData;
	}
	inline unsigned int Capacity() const {
		return alloc_size;
	}
	inline unsigned int Size() const {
		return elemt_size;
	}
	void Copy(const Buffer2D<T>& src,
	          size_t srcBegin, size_t srcEnd,
	          size_t dstBegin) {
		size_t numElmsToCopy = srcEnd - srcBegin;
		if((srcBegin + numElmsToCopy) > alloc_size) {
			/* need to grow */
			ssize_t cap = (srcBegin + numElmsToCopy) - alloc_size;
			ASSERT(cap >= 0);
			Grow(cap + 1000); //sets alloc_size
		}
		memcpy(&data[dstBegin], &src[srcBegin], sizeof(T) * numElmsToCopy);
		ssize_t elemSizeDiff = (ssize_t)(srcBegin + numElmsToCopy) - elemt_size;
		if(elemSizeDiff > 0)
			elemt_size = elemSizeDiff;
	}
private:
	unsigned int elemt_size; //number of objects
	unsigned int alloc_size; //how much is allocated
	T* data;
};
//          x x x
//        x
//0 1 2 3 4 5 6   (7)
#endif
