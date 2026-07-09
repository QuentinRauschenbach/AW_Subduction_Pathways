from dataclasses import dataclass, field

@dataclass
class GateResult:
    data: set
    name: str
    desc: str = ""

    @property
    def length(self):
        return len(self.data)

@dataclass
class GatesSummary:
    name: str
    desc: str = ""
    num_particles: int = 0
    num_aw_particles: int = 0

    # Raw gate crossings
    s1: GateResult = None
    s2: GateResult = None
    s3: GateResult = None
    s4: GateResult = None
    s5: GateResult = None

    # Derived categories
    north:    GateResult = None
    south:    GateResult = None
    yermak:   GateResult = None
    arctic:   GateResult = None
    svalbard: GateResult = None

    # Summary
    all_crossers: GateResult = None
    beached:      GateResult = None
    active:       GateResult = None
    still_active_end_aw: GateResult = None

    @property
    def tracks(self):
        return [self.north, self.south, self.yermak, self.arctic, self.svalbard]
    @property
    def gates(self):
        return [self.s1, self.s2, self.s3, self.s4, self.s5]

    def summary(self):
        gates = [self.north, self.south, self.yermak, self.arctic, 
                 self.svalbard, self.all_crossers, self.beached, self.active]
        width = max(len(g.name) for g in gates if g is not None)
        print(f"\n=== {self.name} ===")
        if self.desc:
            print(f"    {self.desc}")
        if self.num_particles:
            print(f"Total particles: {self.num_particles}")
        if self.num_aw_particles:
            print(f"Particles still active at end of tracking: {self.num_aw_particles}")
        for gate in gates:
            if gate is not None:
                desc_str = f" ({gate.desc})" if gate.desc else ""
                print(f"  {gate.name:{width}} {gate.length:4d}{desc_str}")

        # if neither gate in none print number of not categorizable particles
        if all(g is not None for g in (self.all_crossers, self.north, self.south, self.yermak, self.arctic, self.svalbard)):
            categorized = set().union(*[g.data for g in self.tracks if g is not None])
            uncategorized = self.all_crossers.data - categorized
            print(f"  {'uncategorized':{width}} {len(uncategorized):4d}")

        if all(g is not None for g in (self.beached, self.all_crossers)):
            print(f"  died without crossing: {len(self.beached.data - self.all_crossers.data):3d}")

        if all(g is not None for g in (self.all_crossers, self.still_active_end_aw, self.active)):
            print(f"  still active without crossing (End AW){len(self.active.data - self.all_crossers.data):3d} ({len(self.still_active_end_aw.data - self.all_crossers.data)})")
        
    # Added after saving data files
    @property
    def drifter_to_track(self):
        lookup = {}
        for gate in self.tracks:
            if gate is None:
                continue
            for d_id in gate.data:
                if d_id in lookup:
                    print(f"Warning: Drifter {d_id} found in multiple tracks: {lookup[d_id]} and {gate.name}")
                lookup[d_id] = gate.name
        return lookup