// Copyright (c) 2013-2019 Anton Kozhevnikov, Thomas Schulthess
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that
// the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the
//    following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
//    and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/** \file memory.hpp
 *
 *  \brief Memory management functions and classes.
 */

#ifndef __MEMORY_HPP__
#define __MEMORY_HPP__

#include <list>
#include <iostream>
#include <map>
#include <memory>
#include <cstring>
#include <functional>
#include <algorithm>
#include <array>
#include <complex>
#include <cassert>

#ifdef SIRIUS_USE_MEMORY_POOL
#include <umpire/ResourceManager.hpp>
#include <umpire/Allocator.hpp>
#include <umpire/util/wrap_allocator.hpp>
#include <umpire/strategy/DynamicPoolList.hpp>
#include <umpire/strategy/AlignedAllocator.hpp>
#endif

#include "gpu/acc.hpp"

namespace sddk {

/// Check is the type is a complex number; by default it is not.
template <typename T>
struct is_complex
{
    constexpr static bool value{false};
};

/// Check is the type is a complex number: for std::complex<T> it is true.
template <typename T>
struct is_complex<std::complex<T>>
{
    constexpr static bool value{true};
};

/// Memory types where the code can store data.
/** All memory types can be divided into two (possibly overlapping) groups: accessible by the CPU and accessible by the
 *  device. */
enum class memory_t : unsigned int
{
    /// Nothing.
    none = 0b0000,
    /// Host memory.
    host = 0b0001,
    /// Pinned host memory. This is host memory + extra bit flag.
    host_pinned = 0b0011,
    /// Device memory.
    device = 0b1000,
    /// Managed memory (accessible from both host and device).
    managed = 0b1101,
};

/// Check if this is a valid host memory (memory, accessible by the host).
inline bool
is_host_memory(memory_t mem__)
{
    return static_cast<unsigned int>(mem__) & 0b0001;
}

/// Check if this is a valid device memory (memory, accessible by the device).
inline bool
is_device_memory(memory_t mem__)
{
    return static_cast<unsigned int>(mem__) & 0b1000;
}

/// Get a memory type from a string.
inline memory_t
get_memory_t(std::string name__)
{
    std::transform(name__.begin(), name__.end(), name__.begin(), ::tolower);
    std::map<std::string, memory_t> const m = {{"none", memory_t::none},
                                               {"host", memory_t::host},
                                               {"host_pinned", memory_t::host_pinned},
                                               {"managed", memory_t::managed},
                                               {"device", memory_t::device}};

    if (m.count(name__) == 0) {
        std::stringstream s;
        s << "get_memory_t(): wrong label of the memory_t enumerator: " << name__;
        throw std::runtime_error(s.str());
    }
    return m.at(name__);
}

/// Type of the main processing unit.
/** List the processing units on which the code can run. */
enum class device_t
{
    /// CPU device.
    CPU = 0,

    /// GPU device (with CUDA programming model).
    GPU = 1
};

/// Get type of device by memory type.
inline device_t
get_device_t(memory_t mem__)
{
    switch (mem__) {
        case memory_t::host:
        case memory_t::host_pinned: {
            return device_t::CPU;
        }
        case memory_t::device: {
            return device_t::GPU;
        }
        default: {
            throw std::runtime_error("get_device_t(): wrong memory type");
        }
    }
    return device_t::CPU; // make compiler happy
}

/// Get device type from the string.
inline device_t
get_device_t(std::string name__)
{
    std::transform(name__.begin(), name__.end(), name__.begin(), ::tolower);
    std::map<std::string, device_t> const m = {{"cpu", device_t::CPU}, {"gpu", device_t::GPU}};

    if (m.count(name__) == 0) {
        std::stringstream s;
        s << "get_device_t(): wrong label of the device_t enumerator: " << name__;
        throw std::runtime_error(s.str());
    }
    return m.at(name__);
}

/// Allocate n elements in a specified memory.
/** Allocate a memory block of the memory_t type. Return a nullptr if this memory is not available, otherwise
 *  return a pointer to an allocated block. */
template <typename T>
inline T*
allocate(size_t n__, memory_t M__)
{
    switch (M__) {
        case memory_t::none: {
            return nullptr;
        }
        case memory_t::host: {
            return static_cast<T*>(std::malloc(n__ * sizeof(T)));
        }
        case memory_t::host_pinned: {
#ifdef SIRIUS_GPU
            return acc::allocate_host<T>(n__);
#else
            return nullptr;
#endif
        }
        case memory_t::device: {
#ifdef SIRIUS_GPU
            return acc::allocate<T>(n__);
#else
            return nullptr;
#endif
        }
        default: {
            throw std::runtime_error("allocate(): unknown memory type");
        }
    }
}

/// Deallocate pointer of a given memory type.
inline void
deallocate(void* ptr__, memory_t M__)
{
    switch (M__) {
        case memory_t::none: {
            break;
        }
        case memory_t::host: {
            std::free(ptr__);
            break;
        }
        case memory_t::host_pinned: {
#ifdef SIRIUS_GPU
            acc::deallocate_host(ptr__);
#endif
            break;
        }
        case memory_t::device: {
#ifdef SIRIUS_GPU
            acc::deallocate(ptr__);
#endif
            break;
        }
        default: {
            throw std::runtime_error("deallocate(): unknown memory type");
        }
    }
}

/// Copy between different memory types.
template <typename T, typename F>
inline void
copy(memory_t from_mem__, T const* from_ptr__, memory_t to_mem__, F* to_ptr__, size_t n__)
{
    if (is_host_memory(to_mem__) && is_host_memory(from_mem__)) {
        std::copy(from_ptr__, from_ptr__ + n__, to_ptr__);
        return;
    }
#if defined(SIRIUS_GPU)
    throw std::runtime_error("Copy mixed precision type not supported in device memory");
    return;
#endif
}

template <typename T>
inline void
copy(memory_t from_mem__, T const* from_ptr__, memory_t to_mem__, T* to_ptr__, size_t n__)
{
    if (is_host_memory(to_mem__) && is_host_memory(from_mem__)) {
        std::copy(from_ptr__, from_ptr__ + n__, to_ptr__);
        return;
    }
#if defined(SIRIUS_GPU)
    if (is_device_memory(to_mem__) && is_device_memory(from_mem__)) {
        acc::copy(to_ptr__, from_ptr__, n__);
        return;
    }
    if (is_device_memory(to_mem__) && is_host_memory(from_mem__)) {
        acc::copyin(to_ptr__, from_ptr__, n__);
        return;
    }
    if (is_host_memory(to_mem__) && is_device_memory(from_mem__)) {
        acc::copyout(to_ptr__, from_ptr__, n__);
        return;
    }
#endif
}

/* forward declaration */
class memory_pool;

/// Base class for smart pointer deleters.
class memory_t_deleter_base
{
  protected:
    class memory_t_deleter_base_impl
    {
      public:
        virtual void free(void* ptr__) = 0;
        virtual ~memory_t_deleter_base_impl()
        {
        }
    };
    std::unique_ptr<memory_t_deleter_base_impl> impl_;

  public:
    void operator()(void* ptr__)
    {
        impl_->free(ptr__);
    }
};

/// Deleter for the allocated memory pointer of a given type.
class memory_t_deleter : public memory_t_deleter_base
{
  protected:
    class memory_t_deleter_impl : public memory_t_deleter_base_impl
    {
      protected:
        memory_t M_{memory_t::none};

      public:
        memory_t_deleter_impl(memory_t M__)
            : M_(M__)
        {
        }
        inline void free(void* ptr__)
        {
            if (M_ != memory_t::none) {
                sddk::deallocate(ptr__, M_);
            }
        }
    };

  public:
    explicit memory_t_deleter(memory_t M__)
    {
        impl_ = std::unique_ptr<memory_t_deleter_base_impl>(new memory_t_deleter_impl(M__));
    }
};

/// Deleter for the allocated memory pointer from a given memory pool.
class memory_pool_deleter : public memory_t_deleter_base
{
  protected:
    class memory_pool_deleter_impl : public memory_t_deleter_base_impl
    {
      protected:
        memory_pool* mp_{nullptr};

      public:
        memory_pool_deleter_impl(memory_pool* mp__)
            : mp_(mp__)
        {
        }
        inline void free(void* ptr__);
    };

  public:
    explicit memory_pool_deleter(memory_pool* mp__)
    {
        impl_ = std::unique_ptr<memory_t_deleter_base_impl>(new memory_pool_deleter_impl(mp__));
    }
};

/// Allocate n elements and return a unique pointer.
template <typename T>
inline std::unique_ptr<T, memory_t_deleter_base>
get_unique_ptr(size_t n__, memory_t M__)
{
    return std::unique_ptr<T, memory_t_deleter_base>(allocate<T>(n__, M__), memory_t_deleter(M__));
}

//// Memory pool.
/** This class stores list of allocated memory blocks. Each of the blocks can be divided into subblocks. When subblock
 *  is deallocated it is merged with previous or next free subblock in the memory block. If this was the last subblock
 *  in the block of memory, the (now) free block of memory is merged with the neighbours (if any are available).
 */
class memory_pool
{
  private:
    /// Type of memory that is handeled by this pool.
    memory_t M_;

#ifdef SIRIUS_USE_MEMORY_POOL
    /// handler to umpire allocator_
    umpire::Allocator allocator_;
    /// handler to umpire memory pool
    umpire::Allocator memory_pool_allocator_;
#endif
  public:
    /// Constructor
    memory_pool(memory_t M__, size_t initial_size__ = 0)
        : M_(M__)
    {
        std::string mem_type_;
        // All examples in Umpire use upper case names.

        switch (M__) {
            case memory_t::host: {
                mem_type_ = "HOST";
                break;
            }
            case memory_t::host_pinned: {
                mem_type_ = "PINNED";
                break;
            }
            case memory_t::managed: {
                mem_type_ = "MANAGED";
                break;
            }
            case memory_t::device: {
#ifdef SIRIUS_GPU
                mem_type_ = "DEVICE";
#else
                mem_type_ = "NONE";
                M_ = memory_t::none;
#endif
            }
                break;
            case memory_t::none: {
                mem_type_ = "NONE";
                break;
            }
            default: {
                break;
            }
        }
#ifdef SIRIUS_USE_MEMORY_POOL
        if (M_ != memory_t::none) {
            auto& rm_  = umpire::ResourceManager::getInstance();
            allocator_ = rm_.getAllocator(mem_type_);

            std::transform(mem_type_.begin(), mem_type_.end(), mem_type_.begin(),
                           [](unsigned char c) { return std::tolower(c); });

            if (M_ == memory_t::host) {
                memory_pool_allocator_ =
                    rm_.makeAllocator<umpire::strategy::AlignedAllocator>("aligned_allocator", allocator_, 256);
            } else {
                memory_pool_allocator_ =
                    rm_.makeAllocator<umpire::strategy::DynamicPoolList>(mem_type_ + "_dynamic_pool", allocator_);
            }
        }
#endif
    }

    /// Return a pointer to a memory block for n elements of type T.
    template <typename T>
    T* allocate(size_t num_elements__)
    {
#if defined(SIRIUS_USE_MEMORY_POOL)
        if (M_ == memory_t::none) {
            return nullptr;
        }

        return static_cast<T*>(memory_pool_allocator_.allocate(num_elements__ * sizeof(T)));
#else
        return sddk::allocate<T>(num_elements__, M_);
#endif
    }

    /// Delete a pointer and add its memory back to the pool.
    void free(void* ptr__)
    {
#if defined(SIRIUS_USE_MEMORY_POOL)
        if (M_ == memory_t::none) {
            return;
        }

        memory_pool_allocator_.deallocate(ptr__);
#else
        sddk::deallocate(ptr__, M_);
#endif
    }

    /// Return a unique pointer to the allocated memory.
    template <typename T>
    std::unique_ptr<T, memory_t_deleter_base> get_unique_ptr(size_t n__)
    {
#if defined(SIRIUS_USE_MEMORY_POOL)
        return std::unique_ptr<T, memory_t_deleter_base>(this->allocate<T>(n__), memory_pool_deleter(this));
#else
        return sddk::get_unique_ptr<T>(n__, M_);
#endif
    }

    /// Free all the allocated blocks. umpire does not support this
    /** All pointers and smart pointers, allocated by the pool are invalidated. */
    void reset()
    {
    }

    /// shrink the memory pool and release all memory.
    void clear()
    {
        if (M_ == memory_t::none) {
            return;
        }
#if defined(SIRIUS_USE_MEMORY_POOL)
        memory_pool_allocator_.release();
#endif
    }

    /// Return the type of memory this pool is managing.
    inline memory_t memory_type() const
    {
        return M_;
    }

    /// Return the total capacity of the memory pool.
    size_t total_size() const
    {
#if defined(SIRIUS_USE_MEMORY_POOL)
        if (M_ != memory_t::none) {
            return memory_pool_allocator_.getActualSize();
        }
#endif
        return 0;
    }

    /// Get the total free size of the memory pool.
    size_t free_size() const
    {
#if defined(SIRIUS_USE_MEMORY_POOL)
      if (M_ != memory_t::none) {
          size_t s = memory_pool_allocator_.getActualSize() - memory_pool_allocator_.getCurrentSize();
          return s;
      }
#endif
        return 0;
    }

    /// Get the number of free memory blocks.
    size_t num_blocks() const
    {
#if defined(SIRIUS_USE_MEMORY_POOL)
        if (M_ != memory_t::none) {
            auto dynamic_pool = umpire::util::unwrap_allocator<umpire::strategy::DynamicPoolList>(allocator_);
            return dynamic_pool->getBlocksInPool();
        }
#endif
        return 0;
    }
};
void
memory_pool_deleter::memory_pool_deleter_impl::free(void* ptr__)
{
    mp_->free(ptr__);
}

/// Return a memory pool.
/** A memory pool is created when this function called for the first time. */
sddk::memory_pool& get_memory_pool(sddk::memory_t M__);

#ifdef NDEBUG
#define mdarray_assert(condition__)
#else
#define mdarray_assert(condition__)                                                                                    \
    {                                                                                                                  \
        if (!(condition__)) {                                                                                          \
            std::stringstream _s;                                                                                      \
            _s << "Assertion (" << #condition__ << ") failed "                                                         \
               << "at line " << __LINE__ << " of file " << __FILE__ << std::endl                                       \
               << "array label: " << label_ << std::endl;                                                              \
            for (int i = 0; i < N; i++) {                                                                              \
                _s << "dims[" << i << "].size = " << dims_[i].size() << std::endl;                                     \
            }                                                                                                          \
            throw std::runtime_error(_s.str());                                                                        \
            raise(SIGABRT);                                                                                            \
        }                                                                                                              \
    }
#endif

/// Index descriptor of mdarray.
class mdarray_index_descriptor
{
  public:
    using index_type = int64_t;

  private:
    /// Beginning of index.
    index_type begin_{0};

    /// End of index.
    index_type end_{-1};

    /// Size of index.
    size_t size_{0};

  public:
    /// Constructor of empty descriptor.
    mdarray_index_descriptor()
    {
    }

    /// Constructor for index range [0, size).
    mdarray_index_descriptor(size_t const size__)
        : end_(size__ - 1)
        , size_(size__)
    {
    }

    /// Constructor for index range [begin, end]
    mdarray_index_descriptor(index_type const begin__, index_type const end__)
        : begin_(begin__)
        , end_(end__)
        , size_(end_ - begin_ + 1)
    {
        assert(end_ >= begin_);
    };

    /// Constructor for index range [begin, end]
    mdarray_index_descriptor(std::pair<int, int> const range__)
        : begin_(range__.first)
        , end_(range__.second)
        , size_(end_ - begin_ + 1)
    {
        assert(end_ >= begin_);
    };

    /// Return first index value.
    inline index_type begin() const
    {
        return begin_;
    }

    /// Return last index value.
    inline index_type end() const
    {
        return end_;
    }

    /// Return index size.
    inline size_t size() const
    {
        return size_;
    }

    inline bool check_range(index_type i__) const
    {
#ifdef NDEBUG
        return true;
#else
        if (i__ < begin_ || i__ > end_) {
            std::cout << "index " << i__ << " out of range [" << begin_ << ", " << end_ << "]" << std::endl;
            return false;
        } else {
            return true;
        }
#endif
    }
};

/// Multidimensional array with the column-major (Fortran) order.
/** The implementation supports two memory pointers: one is accessible by CPU and second is accessible by a device.
    The following constructors are implemented:
    \code{.cpp}
    // wrap a host memory pointer and create 2D array 10 x 20.
    mdarray<T, 2>(ptr, 10, 20);

    // wrap a host and device pointers
    mdarray<T, 2>(ptr, ptr_d, 10, 20);

    // wrap a device pointers only
    mdarray<T, 2>(nullptr, ptr_d, 10, 20);

    // create 10 x 20 2D array in main memory
    mdarray<T, 2>(10, 20);

    // create 10 x 20 2D array in device memory
    mdarray<T, 2>(10, 20, memory_t::device);

    // create from the pool memory (pool of any memory type is allowed)
    memory_pool mp(memory_t::host);
    mdarray<T, 2>(mp, 10, 20);
    \endcode
    The pointers can be wrapped only in constructor. Memory allocation can be done by a separate call to .allocate()
    method.
 */
template <typename T, int N>
class mdarray
{
  public:
    using index_type = mdarray_index_descriptor::index_type;

  private:
    /// Optional array label.
    std::string label_;

    /// Unique pointer to the allocated memory.
    std::unique_ptr<T, memory_t_deleter_base> unique_ptr_{nullptr};

    /// Raw pointer.
    T* raw_ptr_{nullptr};
#ifdef SIRIUS_GPU
    /// Unique pointer to the allocated GPU memory.
    std::unique_ptr<T, memory_t_deleter_base> unique_ptr_device_{nullptr};

    /// Raw pointer to GPU memory
    T* raw_ptr_device_{nullptr};
#endif
    /// Array dimensions.
    std::array<mdarray_index_descriptor, N> dims_;

    /// List of offsets to compute the element location by dimension indices.
    std::array<index_type, N> offsets_;

    /// Initialize the offsets used to compute the index of the elements.
    void init_dimensions(std::array<mdarray_index_descriptor, N> const dims__)
    {
        dims_ = dims__;

        offsets_[0] = -dims_[0].begin();
        size_t ld{1};
        for (int i = 1; i < N; i++) {
            ld *= dims_[i - 1].size();
            offsets_[i] = ld;
            offsets_[0] -= ld * dims_[i].begin();
        }
    }

    /// Return linear index in the range [0, size) by the N-dimensional indices (i0, i1, ...)
    template <typename... Args>
    inline index_type idx(Args... args) const
    {
        static_assert(N == sizeof...(args), "wrong number of dimensions");
        std::array<index_type, N> i = {args...};

        for (int j = 0; j < N; j++) {
            // mdarray_assert(dims_[j].check_range(i[j]) && i[j] >= dims_[j].begin() && i[j] <= dims_[j].end());
            mdarray_assert(dims_[j].check_range(i[j]));
        }

        size_t idx = offsets_[0] + i[0];
        for (int j = 1; j < N; j++) {
            idx += i[j] * offsets_[j];
        }
        mdarray_assert(idx >= 0 && idx < size());
        return idx;
    }

    /// Return cosnt pointer to an element at a given index.
    template <bool check_assert = true>
    inline T const* at_idx(memory_t mem__, index_type const idx__) const
    {
        switch (mem__) {
            case memory_t::host:
            case memory_t::host_pinned: {
                if (check_assert) {
                    mdarray_assert(raw_ptr_ != nullptr);
                }
                return &raw_ptr_[idx__];
            }
            case memory_t::device: {
#ifdef SIRIUS_GPU
                if (check_assert) {
                    mdarray_assert(raw_ptr_device_ != nullptr);
                }
                return &raw_ptr_device_[idx__];
#else
                std::printf("error at line %i of file %s: not compiled with GPU support\n", __LINE__, __FILE__);
                throw std::runtime_error("");
#endif
            }
            default: {
                throw std::runtime_error("mdarray::at_idx(): wrong memory type");
            }
        }
        return nullptr; // make compiler happy;
    }

    /// Return pointer to an element at a given index.
    template <bool check_assert = true>
    inline T* at_idx(memory_t mem__, index_type const idx__)
    {
        return const_cast<T*>(static_cast<mdarray<T, N> const&>(*this).at_idx<check_assert>(mem__, idx__));
    }

    // Call constructor on non-trivial data. Complex numbers are treated as trivial.
    inline void call_constructor()
    {
        if (!(std::is_trivial<T>::value || is_complex<T>::value)) {
            for (size_t i = 0; i < size(); i++) {
                new (raw_ptr_ + i) T();
            }
        }
    }

    // Call destructor on non-trivial data. Complex numbers are treated as trivial.
    inline void call_destructor()
    {
        if (!(std::is_trivial<T>::value || is_complex<T>::value)) {
            for (size_t i = 0; i < this->size(); i++) {
                (raw_ptr_ + i)->~T();
            }
        }
    }

    /// Copy constructor is forbidden
    mdarray(mdarray<T, N> const& src) = delete;

    /// Assignment operator is forbidden
    mdarray<T, N>& operator=(mdarray<T, N> const& src) = delete;

  public:
    /// Default constructor.
    mdarray()
    {
    }

    /// Destructor.
    ~mdarray()
    {
        deallocate(memory_t::host);
        deallocate(memory_t::device);
    }

    /// N-dimensional array with index bounds.
    mdarray(std::array<mdarray_index_descriptor, N> const dims__, memory_t memory__ = memory_t::host,
            std::string label__ = "")
        : label_(label__)
    {
        this->init_dimensions(dims__);
        this->allocate(memory__);
    }

    /*
     * 1D array constructors
     *
     */

    /// 1D array with memory allocation.
    mdarray(mdarray_index_descriptor const& d0, memory_t memory__ = memory_t::host, std::string label__ = "")
    {
        static_assert(N == 1, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0});
        this->allocate(memory__);
    }

    /// 1D array with memory pool allocation.
    mdarray(mdarray_index_descriptor const& d0, memory_pool& mp__, std::string label__ = "")
    {
        static_assert(N == 1, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0});
        this->allocate(mp__);
    }

    /// 1D array with host pointer wrapper.
    mdarray(T* ptr__, mdarray_index_descriptor const& d0, std::string label__ = "")
    {
        static_assert(N == 1, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0});
        this->raw_ptr_ = ptr__;
    }

    /// 1D array with host and device pointer wrapper.
    mdarray(T* ptr__, T* ptr_device__, mdarray_index_descriptor const& d0, std::string label__ = "")
    {
        static_assert(N == 1, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0});
        this->raw_ptr_ = ptr__;
#ifdef SIRIUS_GPU
        this->raw_ptr_device_ = ptr_device__;
#endif
    }

    /*
     * 2D array constructors
     *
     */

    /// 2D array with memory allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, memory_t memory__ = memory_t::host,
            std::string label__ = "")
    {
        static_assert(N == 2, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1});
        this->allocate(memory__);
    }

    /// 3D array with memory allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, mdarray_index_descriptor const& d2,
            memory_t memory__ = memory_t::host, std::string label__ = "")
    {
        static_assert(N == 3, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2});
        this->allocate(memory__);
    }

    /// 4D array with memory allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, mdarray_index_descriptor const& d2,
            mdarray_index_descriptor const& d3, memory_t memory__ = memory_t::host, std::string label__ = "")
    {
        static_assert(N == 4, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2, d3});
        this->allocate(memory__);
    }

    /// 5D array with memory allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, mdarray_index_descriptor const& d2,
            mdarray_index_descriptor const& d3, mdarray_index_descriptor const& d4, memory_t memory__ = memory_t::host,
            std::string label__ = "")
    {
        static_assert(N == 5, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2, d3, d4});
        this->allocate(memory__);
    }

    /// 6D array with memory allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, mdarray_index_descriptor const& d2,
            mdarray_index_descriptor const& d3, mdarray_index_descriptor const& d4, mdarray_index_descriptor const& d5,
            memory_t memory__ = memory_t::host, std::string label__ = "")
    {
        static_assert(N == 6, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2, d3, d4, d5});
        this->allocate(memory__);
    }

    /// Wrap a pointer into 2D array.
    mdarray(T* ptr__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, std::string label__ = "")
    {
        static_assert(N == 2, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1});
        this->raw_ptr_ = ptr__;
    }

    mdarray(T* ptr__, T* ptr_device__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1,
            std::string label__ = "")
    {
        static_assert(N == 2, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1});
        this->raw_ptr_ = ptr__;
#ifdef SIRIUS_GPU
        this->raw_ptr_device_ = ptr_device__;
#endif
    }

    /// 2D array with memory pool allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, memory_pool& mp__,
            std::string label__ = "")
    {
        static_assert(N == 2, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1});
        this->allocate(mp__);
    }

    mdarray(T* ptr__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1,
            mdarray_index_descriptor const& d2, std::string label__ = "")
    {
        static_assert(N == 3, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2});
        this->raw_ptr_ = ptr__;
    }

    mdarray(T* ptr__, T* ptr_device__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1,
            mdarray_index_descriptor const& d2, std::string label__ = "")
    {
        static_assert(N == 3, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2});
        this->raw_ptr_ = ptr__;
#ifdef SIRIUS_GPU
        this->raw_ptr_device_ = ptr_device__;
#endif
    }

    /// 3D array with memory pool allocation.
    mdarray(mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1, mdarray_index_descriptor const& d2,
            memory_pool& mp__, std::string label__ = "")
    {
        static_assert(N == 3, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2});
        this->allocate(mp__);
    }

    mdarray(T* ptr__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1,
            mdarray_index_descriptor const& d2, mdarray_index_descriptor const& d3, std::string label__ = "")
    {
        static_assert(N == 4, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2, d3});
        this->raw_ptr_ = ptr__;
    }

    mdarray(T* ptr__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1,
            mdarray_index_descriptor const& d2, mdarray_index_descriptor const& d3, mdarray_index_descriptor const& d4,
            std::string label__ = "")
    {
        static_assert(N == 5, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2, d3, d4});
        this->raw_ptr_ = ptr__;
    }

    mdarray(T* ptr__, mdarray_index_descriptor const& d0, mdarray_index_descriptor const& d1,
            mdarray_index_descriptor const& d2, mdarray_index_descriptor const& d3, mdarray_index_descriptor const& d4,
            mdarray_index_descriptor const& d5, std::string label__ = "")
    {
        static_assert(N == 6, "wrong number of dimensions");

        this->label_ = label__;
        this->init_dimensions({d0, d1, d2, d3, d4, d5});
        this->raw_ptr_ = ptr__;
    }

    /// Move constructor
    mdarray(mdarray<T, N>&& src)
        : label_(src.label_)
        , unique_ptr_(std::move(src.unique_ptr_))
        , raw_ptr_(src.raw_ptr_)
#ifdef SIRIUS_GPU
        , unique_ptr_device_(std::move(src.unique_ptr_device_))
        , raw_ptr_device_(src.raw_ptr_device_)
#endif
    {
        for (int i = 0; i < N; i++) {
            dims_[i]    = src.dims_[i];
            offsets_[i] = src.offsets_[i];
        }
        src.raw_ptr_ = nullptr;
#ifdef SIRIUS_GPU
        src.raw_ptr_device_ = nullptr;
#endif
    }

    /// Move assignment operator
    inline mdarray<T, N>& operator=(mdarray<T, N>&& src)
    {
        if (this != &src) {
            label_       = src.label_;
            unique_ptr_  = std::move(src.unique_ptr_);
            raw_ptr_     = src.raw_ptr_;
            src.raw_ptr_ = nullptr;
#ifdef SIRIUS_GPU
            unique_ptr_device_  = std::move(src.unique_ptr_device_);
            raw_ptr_device_     = src.raw_ptr_device_;
            src.raw_ptr_device_ = nullptr;
#endif
            for (int i = 0; i < N; i++) {
                dims_[i]    = src.dims_[i];
                offsets_[i] = src.offsets_[i];
            }
        }
        return *this;
    }

    /// Allocate memory for array.
    inline mdarray<T, N>& allocate(memory_t memory__)
    {
        /* do nothing for zero-sized array */
        if (!this->size()) {
            return *this;
        }

        /* host allocation */
        if (is_host_memory(memory__)) {
            unique_ptr_ = get_unique_ptr<T>(this->size(), memory__);
            raw_ptr_    = unique_ptr_.get();
            call_constructor();
        }
#ifdef SIRIUS_GPU
        /* device allocation */
        if (is_device_memory(memory__)) {
            unique_ptr_device_ = get_unique_ptr<T>(this->size(), memory__);
            raw_ptr_device_    = unique_ptr_device_.get();
        }
#endif
        return *this;
    }

    /// Allocate memory from the pool.
    inline mdarray<T, N>& allocate(memory_pool& mp__)
    {
        /* do nothing for zero-sized array */
        if (!this->size()) {
            return *this;
        }
        /* host allocation */
        if (is_host_memory(mp__.memory_type())) {
            unique_ptr_ = mp__.get_unique_ptr<T>(this->size());
            raw_ptr_    = unique_ptr_.get();
            call_constructor();
        }
#ifdef SIRIUS_GPU
        /* device allocation */
        if (is_device_memory(mp__.memory_type())) {
            unique_ptr_device_ = mp__.get_unique_ptr<T>(this->size());
            raw_ptr_device_    = unique_ptr_device_.get();
        }
#endif
        return *this;
    }

    /// Deallocate host or device memory.
    inline void deallocate(memory_t memory__)
    {
        if (is_host_memory(memory__)) {
            /* call destructor for non-primitive objects */
            if (unique_ptr_) {
                call_destructor();
            }
            unique_ptr_.reset(nullptr);
            raw_ptr_ = nullptr;
        }
#ifdef SIRIUS_GPU
        if (is_device_memory(memory__)) {
            unique_ptr_device_.reset(nullptr);
            raw_ptr_device_ = nullptr;
        }
#endif
    }

    /// Access operator() for the elements of multidimensional array.
    template <typename... Args>
    inline T const& operator()(Args... args) const
    {
        mdarray_assert(raw_ptr_ != nullptr);
        return raw_ptr_[idx(args...)];
    }

    /// Access operator() for the elements of multidimensional array.
    template <typename... Args>
    inline T& operator()(Args... args)
    {
        return const_cast<T&>(static_cast<mdarray<T, N> const&>(*this)(args...));
    }

    /// Access operator[] for the elements of multidimensional array using a linear index in the range [0, size).
    inline T const& operator[](size_t const idx__) const
    {
        mdarray_assert(idx__ >= 0 && idx__ < size());
        return raw_ptr_[idx__];
    }

    /// Access operator[] for the elements of multidimensional array using a linear index in the range [0, size).
    inline T& operator[](size_t const idx__)
    {
        return const_cast<T&>(static_cast<mdarray<T, N> const&>(*this)[idx__]);
    }

    template <typename... Args>
    inline T const* at(memory_t mem__, Args... args) const
    {
        return at_idx(mem__, idx(args...));
    }

    template <typename... Args>
    inline T* at(memory_t mem__, Args... args)
    {
        return const_cast<T*>(static_cast<mdarray<T, N> const&>(*this).at(mem__, args...));
    }

    /// Return pointer to the beginning of array.
    inline T const* at(memory_t mem__) const
    {
        return at_idx<false>(mem__, 0);
    }

    /// Return pointer to the beginning of array.
    inline T* at(memory_t mem__)
    {
        return const_cast<T*>(static_cast<mdarray<T, N> const&>(*this).at(mem__));
    }

    T* host_data()
    {
        assert(raw_ptr_ != nullptr);
        return raw_ptr_;
    }

    const T* host_data() const
    {
        assert(raw_ptr_ != nullptr);
        return raw_ptr_;
    }

    T* device_data()
    {
#if defined(SIRIUS_GPU)
        assert(raw_ptr_device_ != nullptr);
        return raw_ptr_device_;
#else
        throw std::runtime_error("not compiled with GPU support");
#endif
    }

    const T* device_data() const
    {
#if defined(SIRIUS_GPU)
        assert(raw_ptr_device_ != nullptr);
        return raw_ptr_device_;
#else
        throw std::runtime_error("not compiled with GPU support");
#endif
    }

    /// Return total size (number of elements) of the array.
    inline size_t size() const
    {
        size_t size_{1};

        for (int i = 0; i < N; i++) {
            size_ *= dims_[i].size();
        }

        return size_;
    }

    /// Return size of particular dimension.
    inline size_t size(int i) const
    {
        mdarray_assert(i < N);
        return dims_[i].size();
    }

    /// Return a descriptor of a dimension.
    inline mdarray_index_descriptor dim(int i) const
    {
        mdarray_assert(i < N);
        return dims_[i];
    }

    /// Return leading dimension size.
    inline uint32_t ld() const
    {
        mdarray_assert(dims_[0].size() < size_t(1 << 31));

        return (int32_t)dims_[0].size();
    }

    /// Compute hash of the array
    /** Example: std::printf("hash(h) : %16llX\n", h.hash()); */
    inline uint64_t hash(uint64_t h__ = 5381) const
    {
        for (size_t i = 0; i < size() * sizeof(T); i++) {
            h__ = ((h__ << 5) + h__) + ((unsigned char*)raw_ptr_)[i];
        }

        return h__;
    }

    /// Compute weighted checksum.
    inline T checksum_w(size_t idx0__, size_t size__) const
    {
        T cs{0};
        for (size_t i = 0; i < size__; i++) {
            cs += raw_ptr_[idx0__ + i] * static_cast<double>((i & 0xF) - 8);
        }
        return cs;
    }

    /// Compute checksum.
    inline T checksum(size_t idx0__, size_t size__) const
    {
        T cs{0};
        for (size_t i = 0; i < size__; i++) {
            cs += raw_ptr_[idx0__ + i];
        }
        return cs;
    }

    inline T checksum() const
    {
        return checksum(0, size());
    }

    inline T* begin()
    {
        return this->at(memory_t::host);
    }

    inline T const* begin() const
    {
        return this->at(memory_t::host);
    }

    inline T* end()
    {
        return this->at(memory_t::host) + this->size();
    }

    inline T const* end() const
    {
        return this->at(memory_t::host) + this->size();
    }

    //== template <device_t pu>
    //== inline T checksum() const
    //== {
    //==     switch (pu) {
    //==         case CPU: {
    //==             return checksum();
    //==         }
    //==         case GPU: {
    //==            auto cs = acc::allocate<T>(1);
    //==            acc::zero(cs, 1);
    //==            add_checksum_gpu(raw_ptr_device_, (int)size(), 1, cs);
    //==            T cs1;
    //==            acc::copyout(&cs1, cs, 1);
    //==            acc::deallocate(cs);
    //==            return cs1;
    //==         }
    //==     }
    //== }

    /// Zero n elements starting from idx0.
    inline void zero(memory_t mem__, size_t idx0__, size_t n__)
    {
        mdarray_assert(idx0__ + n__ <= size());
        if (n__ && is_host_memory(mem__)) {
            mdarray_assert(raw_ptr_ != nullptr);
            // std::fill(raw_ptr_ + idx0__, raw_ptr_ + idx0__ + n__, 0);
            std::memset((void*)&raw_ptr_[idx0__], 0, n__ * sizeof(T));
        }
#ifdef SIRIUS_GPU
        if (n__ && on_device() && is_device_memory(mem__)) {
            mdarray_assert(raw_ptr_device_ != nullptr);
            acc::zero(&raw_ptr_device_[idx0__], n__);
        }
#endif
    }

    /// Zero the entire array.
    inline void zero(memory_t mem__ = memory_t::host)
    {
        this->zero(mem__, 0, size());
    }

    /// Copy n elements starting from idx0 from one memory type to another.
    inline void copy_to(memory_t mem__, size_t idx0__, size_t n__, stream_id sid = stream_id(-1))
    {
        if (n__ == 0) {
            return;
        }
#ifdef SIRIUS_GPU
        mdarray_assert(raw_ptr_ != nullptr);
        mdarray_assert(raw_ptr_device_ != nullptr);
        mdarray_assert(idx0__ + n__ <= size());
        if (is_host_memory(mem__) && is_device_memory(mem__)) {
            throw std::runtime_error(
                "mdarray::copy_to(): memory is both host and device, check what to do with this case");
        }
        /* copy to device memory */
        if (is_device_memory(mem__)) {
            if (sid() == -1) {
                /* synchronous (blocking) copy */
                acc::copyin(&raw_ptr_device_[idx0__], &raw_ptr_[idx0__], n__);
            } else {
                /* asynchronous (non-blocking) copy */
                acc::copyin(&raw_ptr_device_[idx0__], &raw_ptr_[idx0__], n__, sid);
            }
        }
        /* copy back from device to host */
        if (is_host_memory(mem__)) {
            if (sid() == -1) {
                /* synchronous (blocking) copy */
                acc::copyout(&raw_ptr_[idx0__], &raw_ptr_device_[idx0__], n__);
            } else {
                /* asynchronous (non-blocking) copy */
                acc::copyout(&raw_ptr_[idx0__], &raw_ptr_device_[idx0__], n__, sid);
            }
        }
#endif
    }

    /// Copy entire array from one memory type to another.
    inline void copy_to(memory_t mem__, stream_id sid = stream_id(-1))
    {
        this->copy_to(mem__, 0, size(), sid);
    }

    /// Check if device pointer is available.
    inline bool on_device() const
    {
#ifdef SIRIUS_GPU
        return (raw_ptr_device_ != nullptr);
#else
        return false;
#endif
    }

    auto label() const
    {
        return label_;
    }

    inline bool on_host() const
    {
        return (raw_ptr_ != nullptr);
    }

    mdarray<T, N>& operator=(std::function<T(void)> f__)
    {
        for (size_t i = 0; i < this->size(); i++) {
            (*this)[i] = f__();
        }
        return *this;
    }

    mdarray<T, N>& operator=(std::function<T(index_type)> f__)
    {
        static_assert(N == 1, "wrong number of dimensions");

        for (index_type i0 = this->dims_[0].begin(); i0 <= this->dims_[0].end(); i0++) {
            (*this)(i0) = f__(i0);
        }
        return *this;
    }

    mdarray<T, N>& operator=(std::function<T(index_type, index_type)> f__)
    {
        static_assert(N == 2, "wrong number of dimensions");

        for (index_type i1 = this->dims_[1].begin(); i1 <= this->dims_[1].end(); i1++) {
            for (index_type i0 = this->dims_[0].begin(); i0 <= this->dims_[0].end(); i0++) {
                (*this)(i0, i1) = f__(i0, i1);
            }
        }
        return *this;
    }
};

// Alias for matrix
template <typename T>
using matrix = mdarray<T, 2>;

/// Serialize to std::ostream
template <typename T, int N>
std::ostream&
operator<<(std::ostream& out, mdarray<T, N> const& v)
{
    if (v.size()) {
        out << v[0];
        for (size_t i = 1; i < v.size(); i++) {
            out << std::string(" ") << v[i];
        }
    }
    return out;
}

/// Copy content of the array to another array of identical size but different precision.
template <typename T, typename F, int N>
inline void
copy(mdarray<F, N> const& src__, mdarray<T, N>& dest__)
{
    if (src__.size() == 0) {
        return;
    }
    for (int i = 0; i < N; i++) {
        if (dest__.dim(i).begin() != src__.dim(i).begin() || dest__.dim(i).end() != src__.dim(i).end()) {
            std::stringstream s;
            s << "error at line " << __LINE__ << " of file " << __FILE__ << " : array dimensions don't match";
            throw std::runtime_error(s.str());
        }
    }
    std::cout << "=== WARNING at line " << __LINE__ << " of file " << __FILE__ << " ===" << std::endl;
    std::cout << "    Copying matrix element with different type, possible loss of data precision" << std::endl;
    std::copy(&src__.at(memory_t::host)[0], &src__.at(memory_t::host)[0] + src__.size(), &dest__.at(memory_t::host)[0]);
}

/// Copy content of the array to another array of identical size.
/** For example:
    \code{.cpp}
    mdarray<double, 2> src(10, 20);
    mdarray<double, 2> dest(10, 20);
    copy(src, dest);
    \endcode
 */
template <typename T, int N>
inline void
copy(mdarray<T, N> const& src__, mdarray<T, N>& dest__)
{
    if (src__.size() == 0) {
        return;
    }
    for (int i = 0; i < N; i++) {
        if (dest__.dim(i).begin() != src__.dim(i).begin() || dest__.dim(i).end() != src__.dim(i).end()) {
            std::stringstream s;
            s << "error at line " << __LINE__ << " of file " << __FILE__ << " : array dimensions don't match";
            throw std::runtime_error(s.str());
        }
    }
    std::copy(&src__[0], &src__[0] + src__.size(), &dest__[0]);
}

/// Copy all memory present on destination.
template <typename T, int N>
void
auto_copy(mdarray<T, N>& dst, const mdarray<T, N>& src)
{

    assert(dst.size() == src.size());
    // TODO: make sure dst and src don't overlap

    if (dst.on_device()) {
        acc::copy(dst.device_data(), src.device_data(), src.size());
    }

    if (dst.on_host()) {
        std::copy(src.host_data(), src.host_data() + dst.size(), dst.host_data());
    }
}

/// Copy memory specified by device from src to dst.
template <typename T, int N>
void
auto_copy(mdarray<T, N>& dst, const mdarray<T, N>& src, device_t device)
{
    // TODO also compare shapes
    if (src.size() == 0) {
        // nothing TODO
        return;
    }

    assert(src.size() == dst.size());
    if (device == device_t::GPU) {
        assert(src.on_device() && dst.on_device());
        acc::copy(dst.device_data(), src.device_data(), dst.size());
    } else if (device == device_t::CPU) {
        assert(src.on_host() && dst.on_host());
        std::copy(src.host_data(), src.host_data() + dst.size(), dst.host_data());
    }
}

template <class numeric_t, std::size_t... Ts>
auto
_empty_like_inner(std::index_sequence<Ts...>& seq, std::size_t (&dims)[sizeof...(Ts)], memory_pool* mempool)
{
    if (mempool == nullptr) {
        return mdarray<numeric_t, sizeof...(Ts)>{dims[Ts]...};
    } else {
        mdarray<numeric_t, sizeof...(Ts)> out{dims[Ts]...};
        out.allocate(*mempool);
        return out;
    }
}

template <typename T, int N>
auto
empty_like(const mdarray<T, N>& src)
{
    auto I = std::make_index_sequence<N>{};
    std::size_t dims[N];
    for (int i = 0; i < N; ++i) {
        dims[i] = src.size(i);
    }
    return _empty_like_inner<T>(I, dims, nullptr);
}

template <typename T, int N>
auto
empty_like(const mdarray<T, N>& src, memory_pool& mempool)
{
    auto I = std::make_index_sequence<N>{};
    std::size_t dims[N];
    for (int i = 0; i < N; ++i) {
        dims[i] = src.size(i);
    }
    return _empty_like_inner<T>(I, dims, &mempool);
}

} // namespace sddk

#endif // __MEMORY_HPP__
