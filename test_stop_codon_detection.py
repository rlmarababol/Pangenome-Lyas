#!/usr/bin/env python3
"""
Self-test for classify_pair() from classify_large_effect_mutations.py,
using synthetic sequences with known expected outcomes. Run from the
same directory as that script so the import resolves.
"""

from classify_large_effect_mutations import classify_pair

# Test 1: stop-gain at codon position 0
# ref codon0=AAA (Lys, sense), codon1=TAA (stop)
# query codon0=TAA (stop), codon1=TAA (stop)
result1 = classify_pair("AAATAA", "TAATAA")
print("Test 1 (expect stop_gain at codon 0):", result1)
assert result1["stop_gain"] == [0], f"FAILED: expected stop_gain=[0], got {result1['stop_gain']}"
assert result1["stop_loss"] == [], f"FAILED: unexpected stop_loss {result1['stop_loss']}"
print("  PASSED\n")

# Test 2: stop-loss at codon position 1
# ref codon0=AAA (sense), codon1=TAA (stop)
# query codon0=AAA (sense), codon1=AAA (sense) -- stop codon replaced by sense
result2 = classify_pair("AAATAA", "AAAAAA")
print("Test 2 (expect stop_loss at codon 1):", result2)
assert result2["stop_loss"] == [1], f"FAILED: expected stop_loss=[1], got {result2['stop_loss']}"
assert result2["stop_gain"] == [], f"FAILED: unexpected stop_gain {result2['stop_gain']}"
print("  PASSED\n")

# Test 3: identical sequences -- no events at all
result3 = classify_pair("AAATAA", "AAATAA")
print("Test 3 (expect no events):", result3)
assert result3["stop_gain"] == [] and result3["stop_loss"] == [], "FAILED: identical seqs should have no stop events"
assert result3["indels"] == [], "FAILED: identical seqs should have no indels"
print("  PASSED\n")

# Test 4: a pair WITH an indel should be skipped entirely for stop-codon
# classification (documented scope), even if it also contains a stop change
result4 = classify_pair("AAATAAGGG", "AAA--AGGG")  # 3bp deletion in query
print("Test 4 (indel present -> stop-codon lists must be empty):", result4)
assert len(result4["indels"]) == 1, f"FAILED: expected 1 indel, got {result4['indels']}"
assert result4["stop_gain"] == [] and result4["stop_loss"] == [], \
    "FAILED: stop-codon classification should be skipped when an indel is present"
print("  PASSED\n")

print("=== All self-tests PASSED. classify_pair() is behaving as designed. ===")
