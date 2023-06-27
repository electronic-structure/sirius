// Copyright (c) 2013-2016 Anton Kozhevnikov, Thomas Schulthess
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

/** \file splindex.hpp
 *
 *  \brief Contains definition of sddk::splindex_base and specializations of sddk::splindex class.
 */

#ifndef __SPLINDEX_HPP__
#define __SPLINDEX_HPP__

#include "strong_type.hpp"
#include "utils/rte.hpp"

/// Basic index type.
template <typename T = int>
struct basic_index_t {
    typedef T value_type;
    typedef T global;
    typedef T local;
};

/// K-point index type.
struct kp_index_t {
    typedef int value_type;
    typedef strong_type<value_type, struct __kp_global_index_tag> global;
    typedef strong_type<value_type, struct __kp_local_index_tag> local;
};

struct atom_index_t {
    typedef int value_type;
    typedef strong_type<value_type, struct __atom_global_index_tag> global;
    typedef strong_type<value_type, struct __atom_local_index_tag> local;
};

struct atom_type_index_t {
    typedef int value_type;
    typedef strong_type<value_type, struct __atom_type_global_index_tag> global;
    typedef strong_type<value_type, struct __atom_type_local_index_tag> local;
};

struct atom_symmetry_class_index_t {
    typedef int value_type;
    typedef strong_type<value_type, struct __atom_symmetry_class_global_index_tag> global;
    typedef strong_type<value_type, struct __atom_symmetry_class_local_index_tag> local;
};

struct paw_atom_index_t {
    typedef int value_type;
    typedef strong_type<value_type, struct __paw_atom_global_index_tag> global;
    typedef strong_type<value_type, struct __paw_atom_local_index_tag> local;
};

/// Number of blocks to which the global index is split.
using n_blocks = strong_type<int, struct __n_blocks_tag>;
/// ID of the block.
/** The id of the block has the range [0, n_blocks) */
using block_id = strong_type<int, struct __block_id_tag>;

namespace sddk {

enum class index_domain_t
{
    /// Global index.
    global,
    /// Local index.
    local
};

/// Base class for split index.
template <typename Index_t = basic_index_t<int>>
class splindex
{
  public:
    using value_type = typename Index_t::value_type;
  protected:
    /// Number of blocks over which the global index is distributed.
    n_blocks n_blocks_{-1};

    /// Index of the block with local fraction of the global index.
    block_id block_id_{-1};

    /// Size (aka length) of the global index.
    value_type size_{-1};

    /// Pair of <local index, block_id> describing the location of a global index element.
    struct location_t
    {
        /// Local index inside a block.
        typename Index_t::local index_local;
        /// Index of the block.
        block_id ib;
        /// Constructor.
        location_t(typename Index_t::local index_local__, block_id ib__)
            : index_local(index_local__)
            , ib(ib__)
        {
        }
    };

  public:

    /// Default constructor.
    splindex()
    {
    }

    /// Constructor.
    /** Check and set index size, number of blocks and block id. */
    splindex(value_type size__, n_blocks n_blocks__, block_id block_id__)
    {
        if (size__ < 0) {
            std::stringstream s;
            s << "wrong size : " << size__;
            throw std::runtime_error(s.str());
        }
        this->size_ = size__;

        if (n_blocks__.get() < 0) {
            std::stringstream s;
            s << "wrong number of blocks : " << n_blocks__.get();
            throw std::runtime_error(s.str());
        }
        this->n_blocks_ = n_blocks__;

        if (block_id__.get() < 0 || block_id__.get() >= n_blocks__.get()) {
            std::stringstream s;
            s << "wrong rank block id : " << block_id__.get();
            throw std::runtime_error(s.str());
        }
        this->block_id_ = block_id__;
    }

    virtual ~splindex()
    {
    }

    /// Return local size of the split index for a given block.
    virtual value_type local_size(block_id block_id__) const = 0;

    /// Return location (block_id and local offset) of the global index.
    virtual location_t location(typename Index_t::global idx__) const = 0;

    /// Return global index by block id and local index.
    virtual typename Index_t::global global_index(typename Index_t::local idxloc__, block_id block_id__) const = 0;

    /// Return local size for the current block.
    value_type local_size() const
    {
        return this->local_size(this->block_id_);
    }

    /// Return global index of an element by local index and block id.
    inline auto global_index(typename splindex<Index_t>::location_t loc__) const
    {
        return this->global_index(loc__.index_local, loc__.ib);
    }

    inline auto global_index(typename Index_t::local idxloc__) const
    {
        return this->global_index(idxloc__, this->block_id_);
    }

    /// Return total length of the index (global number of elements).
    inline auto size() const noexcept
    {
        return size_;
    }

    /// Compute size of the block from global index size and number of blocks.
    static inline auto block_size(value_type size__, n_blocks n_blocks__)
    {
        return size__ / n_blocks__ + std::min(value_type(1), size__ % n_blocks__);
    }
};

template <typename Index_t>
class splindex_iterator_t : public std::iterator<std::random_access_iterator_tag, Index_t>
{
  private:
    splindex<Index_t> const& idx_;
  public:
    using difference_type = typename std::iterator<std::random_access_iterator_tag, Index_t>::difference_type;
    typename Index_t::local li;
    typename Index_t::global i;

    splindex_iterator_t<Index_t>& operator=(splindex_iterator_t<Index_t> const& lhs_) = default;
    //{
    //    this->li = lhs_.li;
    //    this->i = lhs_.i;
    //    return *this;
    //}

    splindex_iterator_t(splindex<Index_t> const& idx__)
        : idx_{idx__}
        , li{0}
        , i{0}
    {
    }
    inline bool operator!=(splindex_iterator_t<Index_t> const& rhs__)
    {
        return this->li != rhs__.li;
    }
    inline splindex_iterator_t<Index_t>& operator++()
    {
        this->li++;
        return *this;
    }
    inline splindex_iterator_t<Index_t>& operator++(int)
    {
        splindex_iterator_t<Index_t> tmp(this->idx());
        this->li++;
        return tmp;
    }
    inline splindex_iterator_t<Index_t> const& operator*()
    {
        this->i = idx_.global_index(this->li);
        return *this;
    }
    inline difference_type operator-(splindex_iterator_t<Index_t> const& rhs__) const
    {
        return li - rhs__.li;
    }
    inline splindex_iterator_t<Index_t>& operator+=(difference_type rhs__)
    {
        li += rhs__;
        return *this;
    }
};

template <typename Index_t = basic_index_t<int>>
class splindex_block : public splindex<Index_t>
{
  public:
    using value_type = typename splindex<Index_t>::value_type;
  private:
    /// Local index size of a given block.
    value_type block_size_;
  public:
    splindex_block()
    {
    }
    /// Constructor.
    splindex_block(value_type size__, n_blocks n_blocks__, block_id block_id__)
        : splindex<Index_t>(size__, n_blocks__, block_id__)
    {
        this->block_size_ = this->block_size(size__, n_blocks__);
    }

    using splindex<Index_t>::local_size;

    /// Return local size of the split index for a given block.
    inline value_type local_size(block_id block_id__) const
    {
        RTE_ASSERT(block_id__ >= 0 && block_id__ < this->n_blocks_);

        if (this->size_ == 0) {
            return 0;
        }

        auto n = static_cast<int>(this->size_ / block_size_);
        if (block_id__ < n) {
            return block_size_;
        } else {
            return std::max(0, this->size_ - block_id__ * block_size_);
        }
    }

    /// Return "local index, rank" pair for a global index.
    inline typename splindex<Index_t>::location_t location(typename Index_t::global idx__) const
    {
        RTE_ASSERT(idx__ < this->size_);

        auto ib = static_cast<int>(idx__ / this->block_size_);
        value_type idxloc = idx__ - ib * this->block_size_;

        return typename splindex<Index_t>::location_t(typename Index_t::local(idxloc), block_id(ib));
    }

    using splindex<Index_t>::global_index;

    /// Return global index of an element by local index and block id.
    inline typename Index_t::global global_index(typename Index_t::local idxloc__, block_id block_id__) const
    {
        RTE_ASSERT(block_id__ >= 0 && block_id__ < this->n_blocks_);

        if (this->local_size(block_id__) == 0) {
            return typename Index_t::global(-1);
        }

        RTE_ASSERT(idxloc__ < local_size(block_id__));

        return typename Index_t::global(this->block_size_ * block_id__ + idxloc__);
    }

    inline auto global_offset() const
    {
        return this->global_index(typename Index_t::local(0), this->block_id_);
    }

    inline auto global_offset(block_id iblock__) const
    {
        return this->global_index(typename Index_t::local(0), iblock__);
    }

    inline auto counts() const
    {
        std::vector<value_type> v(this->n_blocks_);
        for (int i = 0; i < this->n_blocks_; i++) {
            v[i] = local_size(block_id(i));
        }
        return v;
    }
};

template <typename Index_t = basic_index_t<int>>
class splindex_block_cyclic : public splindex<Index_t>
{
  public:
    using value_type = typename splindex<Index_t>::value_type;
  private:
    /// Cyclic block size.
    value_type block_size_;
  public:
    splindex_block_cyclic()
    {
    }
    /// Constructor.
    splindex_block_cyclic(value_type size__, n_blocks n_blocks__, block_id block_id__, value_type block_size__)
        : splindex<Index_t>(size__, n_blocks__, block_id__)
        , block_size_{block_size__}
    {
    }

    using splindex<Index_t>::local_size;

    /// Return local size of the split index for a given block.
    inline value_type local_size(block_id block_id__) const
    {
        RTE_ASSERT(block_id__ >= 0 && block_id__ < this->n_blocks_);

        if (this->size_ == 0) {
            return 0;
        }
        /* number of full blocks */
        auto num_blocks = this->size() / this->block_size_;

        auto n = (num_blocks / this->n_blocks_) * this->block_size_;
        auto rank_offs = static_cast<int>(num_blocks % this->n_blocks_);

        if (block_id__ < rank_offs) {
            n += this->block_size_;
        } else if (block_id__ == rank_offs) {
            n += this->size_ % this->block_size_;
        }
        return n;
    }

    /// Return "local index, rank" pair for a global index.
    inline typename splindex<Index_t>::location_t location(typename Index_t::global idx__) const
    {
        RTE_ASSERT(idx__ < this->size_);

        /* number of full blocks */
        auto num_blocks = idx__ / this->block_size_;

        /* local index */
        value_type idxloc = (num_blocks / this->n_blocks_) * block_size_ + idx__ % this->block_size_;

        /* corresponding rank */
        auto ib = static_cast<int>(num_blocks % this->n_blocks_);

        return typename splindex<Index_t>::location_t(typename Index_t::local(idxloc), block_id(ib));
    }

    using splindex<Index_t>::global_index;

    /// Return global index of an element by local index and block id.
    inline typename Index_t::global global_index(typename Index_t::local idxloc__, block_id block_id__) const
    {
        RTE_ASSERT(block_id__ >= 0 && block_id__ < this->n_blocks_);
        RTE_ASSERT(idxloc__ < local_size(block_id__));

        auto nb = idxloc__ / this->block_size_;

        return typename Index_t::global((nb * this->n_blocks_ + block_id__) * this->block_size_ +
                idxloc__ % this->block_size_);
    }

};

/// Externally defined block distribution.
template <typename Index_t = basic_index_t<int>>
class splindex_chunk : public splindex<Index_t>
{
  public:
    using value_type = typename splindex<Index_t>::value_type;
  private:
    std::vector<std::vector<value_type>> global_index_;
    std::vector<typename splindex<Index_t>::location_t> locations_;

  public:
    /// Default constructor.
    splindex_chunk()
    {
    }

    /// Constructor with specific partitioning.
    splindex_chunk(value_type size__, n_blocks n_blocks__, block_id block_id__, std::vector<value_type> const counts__)
        : splindex<Index_t>(size__, n_blocks__, block_id__)
    {
        for (int r = 0; r < n_blocks__.get(); r++) {
            global_index_.push_back(std::vector<value_type>());
            for (value_type i = 0; i < counts__[r]; i++) {
                global_index_.back().push_back(static_cast<value_type>(locations_.size()));
                locations_.push_back(typename splindex<Index_t>::location_t(typename Index_t::local(i), block_id(r)));
            }
        }

        RTE_ASSERT(static_cast<value_type>(locations_.size()) == this->size());
    }

    using splindex<Index_t>::local_size;

    inline value_type local_size(block_id block_id__) const
    {
        RTE_ASSERT(block_id__ >= 0 && block_id__ < this->n_blocks_);

        if (this->size_ == 0) {
            return 0;
        }
        return static_cast<value_type>(global_index_[block_id__].size());
    }

    inline typename splindex<Index_t>::location_t location(typename Index_t::global idx__) const
    {
        return locations_[idx__];
    }

    using splindex<Index_t>::global_index;

    inline typename Index_t::global global_index(typename Index_t::local idxloc__, block_id block_id__) const
    {
        RTE_ASSERT(block_id__ >= 0 && block_id__ < this->n_blocks_);
        RTE_ASSERT(idxloc__ < local_size(block_id__));

        return typename Index_t::global(global_index_[block_id__][idxloc__]);
    }

    inline auto global_offset() const
    {
        return this->global_index(typename Index_t::local(0), this->block_id_);
    }
};

template <typename Index_t>
auto begin_global(splindex<Index_t> const& a__)
{
    return typename Index_t::global(0);
}

template <typename Index_t>
auto end_global(splindex<Index_t> const& a__)
{
    return typename Index_t::global(a__.size());
}

template <typename Index_t>
auto begin(splindex<Index_t> const& a__)
{
    splindex_iterator_t<Index_t> it(a__);
    it.li = typename Index_t::local(0);
    return it;
}

template <typename Index_t>
auto end(splindex<Index_t> const& a__)
{
    splindex_iterator_t<Index_t> it(a__);
    it.li = typename Index_t::local(a__.local_size());
    return it;
}

} // namespace sddk

#endif // __SPLINDEX_HPP__
