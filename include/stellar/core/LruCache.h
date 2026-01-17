#pragma once

#include <cstddef>
#include <functional>
#include <list>
#include <unordered_map>
#include <utility>

namespace stellar::core {

// A small, general-purpose LRU cache.
//
// - O(1) find/insert/erase.
// - Most-recently-used entries are kept at the front of an internal list.
// - Capacity == 0 means "unbounded" (no automatic eviction).
//
// Intended for caches that sit above deterministic/procedural generators,
// expensive parsers, etc.
template <typename Key,
          typename Value,
          typename Hash = std::hash<Key>,
          typename Eq = std::equal_to<Key>>
class LruCache {
public:
  explicit LruCache(std::size_t capacity = 0) : capacity_(capacity) {}

  std::size_t capacity() const { return capacity_; }
  std::size_t size() const { return map_.size(); }
  bool empty() const { return map_.empty(); }

  void reserve(std::size_t n) { map_.reserve(n); }

  // Set maximum number of entries. capacity == 0 => unbounded.
  void setCapacity(std::size_t cap) {
    capacity_ = cap;
    pruneToCapacity();
  }

  // Same as setCapacity(), but calls onEvict(key, value) for each removed entry.
  // Returns number of evicted entries.
  template <typename EvictFn>
  std::size_t setCapacity(std::size_t cap, EvictFn&& onEvict) {
    capacity_ = cap;
    return pruneToCapacity(std::forward<EvictFn>(onEvict));
  }

  // Clear all entries.
  void clear() {
    map_.clear();
    order_.clear();
  }

  // Find an entry and mark it as most-recently used.
  Value* find(const Key& key) {
    auto it = map_.find(key);
    if (it == map_.end()) return nullptr;
    touch(it);
    return &it->second.value;
  }

  // Peek without updating recency.
  const Value* peek(const Key& key) const {
    auto it = map_.find(key);
    if (it == map_.end()) return nullptr;
    return &it->second.value;
  }

  bool contains(const Key& key) const { return map_.find(key) != map_.end(); }

  // Erase an entry. Returns true if the key existed.
  bool erase(const Key& key) {
    auto it = map_.find(key);
    if (it == map_.end()) return false;
    order_.erase(it->second.it);
    map_.erase(it);
    return true;
  }

  // Get an entry if present; otherwise create it with createFn().
  //
  // createFn must be callable with no args and return a Value.
  template <typename CreateFn>
  Value& getOrCreate(const Key& key, CreateFn&& createFn, bool* outCreated = nullptr) {
    return getOrCreate(key, std::forward<CreateFn>(createFn), outCreated, [](const Key&, Value&) {});
  }

  // Same as getOrCreate(), but calls onEvict(key, value) for each removed entry.
  template <typename CreateFn, typename EvictFn>
  Value& getOrCreate(const Key& key,
                     CreateFn&& createFn,
                     bool* outCreated,
                     EvictFn&& onEvict) {
    auto it = map_.find(key);
    if (it != map_.end()) {
      touch(it);
      if (outCreated) *outCreated = false;
      return it->second.value;
    }

    // Insert at MRU position.
    order_.push_front(key);
    const auto lit = order_.begin();

    Entry e{std::forward<CreateFn>(createFn)(), lit};
    auto [mit, ok] = map_.emplace(key, std::move(e));
    if (!ok) {
      // Extremely unlikely (another emplace won a race in a multithreaded
      // scenario). Fall back to existing entry.
      order_.pop_front();
      touch(mit);
      if (outCreated) *outCreated = false;
      return mit->second.value;
    }

    if (outCreated) *outCreated = true;

    pruneToCapacity(std::forward<EvictFn>(onEvict));
    return mit->second.value;
  }

  // Insert a value if the key is absent; otherwise return the existing entry.
  //
  // This is useful when value construction is expensive and you want to perform it
  // *outside* the cache lock (e.g., parallel generation), then insert the ready
  // value and still preserve LRU behavior.
  template <typename EvictFn>
  Value& getOrInsert(const Key& key,
                     Value value,
                     bool* outInserted,
                     EvictFn&& onEvict) {
    auto it = map_.find(key);
    if (it != map_.end()) {
      touch(it);
      if (outInserted) *outInserted = false;
      return it->second.value;
    }

    order_.push_front(key);
    const auto lit = order_.begin();

    Entry e{std::move(value), lit};
    auto [mit, ok] = map_.emplace(key, std::move(e));
    if (!ok) {
      order_.pop_front();
      touch(mit);
      if (outInserted) *outInserted = false;
      return mit->second.value;
    }

    if (outInserted) *outInserted = true;
    pruneToCapacity(std::forward<EvictFn>(onEvict));
    return mit->second.value;
  }

  Value& getOrInsert(const Key& key, Value value, bool* outInserted = nullptr) {
    return getOrInsert(key, std::move(value), outInserted, [](const Key&, Value&) {});
  }

  // Insert or replace an entry. Marks the entry as MRU.
  // Returns true if inserted new, false if replaced.
  template <typename EvictFn>
  bool insertOrAssign(const Key& key, Value value, EvictFn&& onEvict) {
    auto it = map_.find(key);
    if (it != map_.end()) {
      it->second.value = std::move(value);
      touch(it);
      return false;
    }

    order_.push_front(key);
    const auto lit = order_.begin();

    Entry e{std::move(value), lit};
    map_.emplace(key, std::move(e));

    pruneToCapacity(std::forward<EvictFn>(onEvict));
    return true;
  }

  bool insertOrAssign(const Key& key, Value value) {
    return insertOrAssign(key, std::move(value), [](const Key&, Value&) {});
  }

  // Evict least-recently-used entries until size() <= capacity().
  // Returns number of evicted entries.
  std::size_t pruneToCapacity() {
    return pruneToCapacity([](const Key&, Value&) {});
  }

  // Same as pruneToCapacity(), but calls onEvict(key, value) for each removed entry.
  template <typename EvictFn>
  std::size_t pruneToCapacity(EvictFn&& onEvict) {
    if (capacity_ == 0) return 0; // unbounded

    std::size_t evicted = 0;
    while (map_.size() > capacity_) {
      auto last = std::prev(order_.end());
      const Key& k = *last;

      auto it = map_.find(k);
      if (it != map_.end()) {
        std::forward<EvictFn>(onEvict)(it->first, it->second.value);
        map_.erase(it);
      }

      order_.erase(last);
      ++evicted;
    }
    return evicted;
  }

private:
  using List = std::list<Key>;
  using ListIt = typename List::iterator;

  struct Entry {
    Value value;
    ListIt it;
  };

  using Map = std::unordered_map<Key, Entry, Hash, Eq>;
  using MapIt = typename Map::iterator;

  void touch(MapIt it) {
    // Move this key to the front (MRU). Splice is O(1) and does not invalidate iterators.
    order_.splice(order_.begin(), order_, it->second.it);
  }

  std::size_t capacity_{0};
  List order_{};
  Map map_{};
};

} // namespace stellar::core
