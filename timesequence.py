from typing import List, Tuple
from itertools import combinations

def generate_dispatch_schedules(
    slots1: List[float],
    slots2: List[float],
    total_serviceable: int,
    minQ: int,
    maxQ: int,
    Xmin: int = 4,
    Xmax: int = 10
) -> List[Tuple[List[float], List[int]]]:
    merged_slots = sorted(slots1 + slots2)
    N = len(merged_slots)
    results: List[Tuple[List[float], List[int]]] = []

    # never try to pick more slots than you have
    Xmax = min(Xmax, N)

    def backtrack_quantities(
        idx: int,
        x: int,
        remaining: int,
        current_Qs: List[int],
        current_times: List[float]
    ):
        if idx == x:
            if remaining == 0:
                # Found one
                results.append((current_times.copy(), current_Qs.copy()))
            return

        slots_left = x - idx
        # we will pick one q now, then slots_left-1 remain
        for q in range(minQ, maxQ + 1):
            rem_after = remaining - q
            min_possible = (slots_left - 1) * minQ
            max_possible = (slots_left - 1) * maxQ
            # prune
            if rem_after < min_possible or rem_after > max_possible:
                continue
            current_Qs.append(q)
            backtrack_quantities(idx + 1, x, rem_after, current_Qs, current_times)
            current_Qs.pop()

    def backtrack_times(
        idx: int,
        x: int,
        start_i: int,
        current_times: List[float]
    ):
        if idx == x:
            # got x times, now assign quantities
            backtrack_quantities(0, x, total_serviceable, [], current_times)
            return

        # we need to leave room for (x-idx) picks total
        for i in range(start_i, N - (x - idx) + 1):
            current_times.append(merged_slots[i])
            backtrack_times(idx + 1, x, i + 1, current_times)
            current_times.pop()

    for x in range(Xmin, Xmax + 1):
        print(f"--- Trying x = {x} dispatches ---")
        backtrack_times(0, x, 0, [])

    return results


if __name__ == "__main__":
    # your two windows, discretized however you like:
    slots1 = list(range(4, 10))    # e.g. T1=1, T3=3
    slots2 = list(range(17, 21))   # e.g. T5=5,  T =10

    total_serviceable = 200
    minQ = 10
    maxQ = 15

    schedules = generate_dispatch_schedules(
        slots1, slots2,
        total_serviceable,
        minQ, maxQ,
        Xmin=4, Xmax=10
    )
    print(f"\nFound {len(schedules)} schedules\n")

    # print the first 10 to check
    for times, qs in schedules[:10]:
        print(f"Times: {times}  Qs: {qs}")