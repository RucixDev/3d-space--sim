#include "test_harness.h"

#include "stellar/core/LruCache.h"

#include <string>

int test_lru_cache() {
  int failures = 0;

  {
    stellar::core::LruCache<int, std::string> c(2);
    bool created = false;

    auto& v1 = c.getOrCreate(1, [] { return std::string("one"); }, &created);
    CHECK(created);
    CHECK(v1 == "one");

    c.getOrCreate(2, [] { return std::string("two"); }, &created);
    CHECK(created);
    CHECK(c.size() == 2);

    // Touch key=1 so key=2 becomes LRU.
    CHECK(c.find(1) != nullptr);

    c.getOrCreate(3, [] { return std::string("three"); }, &created);
    CHECK(created);
    CHECK(c.size() == 2);

    // key=2 should have been evicted.
    CHECK(c.find(2) == nullptr);
    CHECK(c.find(1) != nullptr);
    CHECK(c.find(3) != nullptr);
  }

  {
    // getOrInsert() inserts a value that was created outside the cache.
    stellar::core::LruCache<int, std::string> c(2);
    bool inserted = false;

    auto& v = c.getOrInsert(1, std::string("one"), &inserted);
    CHECK(inserted);
    CHECK(v == "one");

    // Calling again should return the existing entry and not replace the value.
    auto& v2 = c.getOrInsert(1, std::string("uno"), &inserted);
    CHECK(!inserted);
    CHECK(v2 == "one");

    // insertOrAssign() replaces the value (and returns false because it was not a new insert).
    const bool insNew = c.insertOrAssign(1, std::string("uno"));
    CHECK(!insNew);
    CHECK(c.peek(1) && *c.peek(1) == "uno");
  }

  {
    // Evict callback works for getOrInsert().
    stellar::core::LruCache<int, int> c(1);
    bool inserted = false;

    CHECK(c.getOrInsert(1, 11, &inserted) == 11);
    CHECK(inserted);

    int evictedKey = -1;
    int evictedValue = -1;
    c.getOrInsert(2, 22, &inserted, [&](const int& k, int& v) {
      evictedKey = k;
      evictedValue = v;
    });
    CHECK(inserted);
    CHECK(evictedKey == 1);
    CHECK(evictedValue == 11);
  }

  {
    // capacity=0 => unbounded.
    stellar::core::LruCache<int, int> c(0);
    for (int i = 0; i < 1000; ++i) {
      c.getOrCreate(i, [i] { return i * 2; });
    }
    CHECK(c.size() == 1000);
    CHECK(c.peek(123) && *c.peek(123) == 246);
  }

  return failures;
}
