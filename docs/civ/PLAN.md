# `fms-civ` plan

The canonical design of `fms-civ`: what it is for and the decisions that shape
it. Settled on 2026-11-12; nothing is built yet. The family it belongs to is in
[../PLAN.md](../PLAN.md) §1, and the file format it reads and writes is in
[../FMSMAP.md](../FMSMAP.md).

`fms-civ` is where you put civilisations into a world and let them loose: draw
political borders, give each faction a few traits, and watch them grow, shrink,
fight and split on the land the other modules made. It is a tool for roleplay, a
very limited god game, and not a strategy game or a model of real history.

## 1. Goals

- **Imagine your civilisations in your world**: the pitch is roleplay, not
  accurate history.
- **Simple by design**: no technology, no strategy layer, nothing to win.
- **The land matters**: factions on poor land shrink or die out, and shortage of
  food leads to war. This is what the climate and geology simulations are for.
- **Factions, cultures and religions**: three maps that change over time.
- **A limited god game**: disasters happen naturally and can also be triggered
  by hand.
- **Stories out**: a scrubbable animated map, a written chronicle, graphs, and
  layers back into `.fmsmap`.

## 2. Key decisions

| Area | Decision |
|---|---|
| Simulated unit | Provinces, each with a population, a food capacity, culture and religion shares, and a main town |
| Provinces | Generated along rivers, ridges and coasts to a target count; borders you draw become province edges; provinces can be merged, split and renamed |
| Resource | Food capacity only |
| Faction traits | A few sliders the simulation reads, plus name, colour and free lore text it ignores |
| Culture and religion | Their own layers, independent of political borders, held as shares per province; one effect: breakaways and war are likelier between unlike people |
| Conflict | Wars as dated, named events; no armies |
| Unclaimed land | Empty until settled |
| Sea | Nobody owns it; fishing on coasts and short crossings |
| Disasters | Eruptions and earthquakes, storms and floods, droughts and cold years; natural or triggered by hand |
| Time | One step per year; typical runs of 200–2,000 years |
| Interaction | A run takes no input while playing; at pause you can edit anything or trigger a disaster |
| Timeline | Every year recorded and scrubbable; changing an earlier year discards what followed |
| Names | Generated per culture from the seed; anything can be renamed |
| The world | Read-only; civ writes only its own layers |
| Outputs | Animated map with timeline, chronicle, per-faction graphs, layers in `.fmsmap` with name tables |
| Application | Its own GUI program on the shared render crate, with the simulation as a headless library |
| Determinism | Seed, setup and the log of edits and disasters reproduce a history exactly |
| Validation | Behaviour tests, plus one look at Earth for tuning |
| First milestone | Headless core |

## 3. The world as civ sees it

### 3.1 Provinces

A province is a set of cells of the terrain grid. civ generates them along
rivers, ridges and coasts, smaller where the land is rich, to a target count
that is a setting (a few thousand by default on an Earth-sized planet). Any
border you draw is also a province edge, so drawn borders are exact. Provinces
can be merged, split and renamed. The sea is not cut into provinces. The
simulation runs per province; maps are drawn per cell.

### 3.2 Food capacity

The only resource: how many people a province's land can feed. It is computed
from whatever layers the file has.

| Input | Source |
|---|---|
| Temperature and rainfall | `fms-climate` results; else the biome → typical-climate table, with its warning ([../FMSMAP.md](../FMSMAP.md)); else a crude guess from latitude and altitude, with a stronger warning |
| Rivers and lakes | The file's river and lake layers |
| Coast | Fishing in coastal provinces |
| Slope and altitude | From `terrain.elevation` |
| Soil | A bonus from `fms-geo`'s sediment thickness and volcanic rock, when present |

The capacity a faction gets from a province also depends on its way of life
(§4.1). The world is read-only: farming, clearing and irrigation are not
modelled, and civ never changes a layer it did not write.

### 3.3 Towns

Each province has one named town on its best cell (river, coast, richest land).
Its size follows the province's population. Each faction's capital is one of its
towns.

## 4. Factions

### 4.1 Traits

A faction has a name, a colour, free lore text the simulation ignores, a
capital, and a few sliders it reads. The intended sliders (exact list fixed when
built):

| Trait | Effect |
|---|---|
| Expansionism | How readily spare people settle new provinces |
| Militancy | How readily it goes to war, and part of its strength |
| Cohesion | How well it holds distant or unlike provinces |
| Way of life | Farming, herding or seafaring: which land feeds it best, and how far it reaches over water |

### 4.2 Population and settling

Each province's population grows toward its food capacity. Shortage means
famine and decline, so a faction on land that cannot feed it shrinks and can die
out. A faction with people to spare settles neighbouring empty provinces, or
ones across a short stretch of water. Land outside any faction is empty until
settled. Growth rates and the settling rule are fixed when built.

### 4.3 War

Wars are events with a start, an end and a name. A faction short of food
declares war on a neighbour, more readily if it is militant or the neighbour is
unlike it in culture or religion. While the war lasts, border provinces change
hands toward the stronger side. Strength comes from population and militancy;
distance from the capital and rough terrain favour the defender. Then peace. A
faction with no provinces left is gone. There are no armies, tactics or
alliances. The exact odds and durations are fixed when built.

### 4.4 Breakaways

Provinces whose people are unlike their rulers in culture or religion, or that
lie far from the capital, may break away as a new faction, with traits inherited
from the old one and a generated name. Low cohesion makes this likelier. A
faction can shatter altogether; provinces left without a faction stay populated
and stateless until another faction absorbs them.

### 4.5 The sea

Nobody owns the sea. Coastal provinces get food from fishing, and factions
settle or attack across short crossings; a seafaring way of life lengthens the
reach. There are no ships or navies.

## 5. Culture and religion

Two layers of their own, independent of political borders. A province holds
shares (say 70% one culture, 30% another); the maps show the majority. Shares
drift toward neighbouring provinces and toward the rulers, culture slowly and
religion faster. Their single effect on the rest of the simulation is the one in
§4.3 and §4.4. Only the cultures and religions you define exist; none arise by
themselves in v1. The drift rates are fixed when built.

## 6. Disasters

| Disaster | Natural rate and place from |
|---|---|
| Eruptions and earthquakes | `fms-geo`'s volcanism and seismic-hazard layers |
| Storms and floods | `fms-climate`'s storm and flood statistics |
| Droughts and cold years | `fms-climate`'s year-to-year swings in rain and temperature |

Natural disasters strike by seeded chance at those rates; when the layers are
missing civ uses defaults and warns. You can also trigger any of them by hand,
where you choose. A disaster kills people and cuts food capacity in the
provinces it hits for a while, and goes into the chronicle. `fms-climate` has to
export the statistics named here ([../FMSMAP.md](../FMSMAP.md)).

## 7. Running and editing

### 7.1 Setup

Draw borders, define factions and their traits, paint cultures and religions,
set the seed, province count and run length. Starting population follows from
food capacity unless set by hand.

### 7.2 A run

One step per year; typical runs cover 200–2,000 years. While playing, a run
takes no input. At pause you can do anything setup allows (redraw a border,
change a trait, add or remove a faction, repaint culture) or trigger a disaster.
Every such change is logged.

### 7.3 Timeline

Every year is recorded, so the map can be scrubbed back and forth. Changing
something in an earlier year discards everything after it and the run continues
from there. A project has one timeline; keeping the old future means saving a
copy first.

### 7.4 Names

Each culture gets a simple sound pattern from the seed, so its names sound
alike. Provinces, towns, breakaway factions, wars and disasters are named from
it. Anything can be renamed, and a name you typed is never overwritten.

### 7.5 Determinism

The seed, the setup and the log of edits and triggered disasters reproduce a
history exactly; a new seed gives an alternate history. The cross-platform
guarantee is the family's ([../PLAN.md](../PLAN.md) §2).

## 8. Outputs

- **Animated map** of factions, cultures, religions and population, with a
  timeline to scrub.
- **Chronicle**: a dated list of events (settlings, wars, provinces and capitals
  lost, breakaways, collapses, famines, disasters), exportable as text.
- **Graphs**: population and territory per faction over time.
- **`.fmsmap` layers** at a chosen year (§9).

## 9. Files

### 9.1 Layers written

Names and encodings are fixed when each is built.

| Layer | Notes |
|---|---|
| Province | Id per cell |
| Faction | Id per cell; none for unclaimed or stateless land |
| Culture | Majority per cell |
| Religion | Majority per cell |
| Population density | People per area |
| Towns | The cell of each province's town |

Each id layer comes with a small table giving every id its name and colour, a
format addition ([../FMSMAP.md](../FMSMAP.md)). Label placement and styling stay
in `fms-editor`, which can draw its political maps from these layers. A later
civ run can start from them: civ's own layers are treated as initial. Every
layer civ does not write passes through unchanged.

### 9.2 Project file

civ's own file holds what other modules never read: setup, seed and settings,
the recorded timeline, the chronicle, and the log of edits and triggered
disasters.

## 10. Application

A separate `fms-civ` program built on the render crate it shares with
`fms-viewer` and `fms-editor`, with the crates in the family workspace. The
simulation is a library with no GUI, with a headless runner for tests and
scripted setups. The GUI adds the drawing tools, trait panels, the timeline, the
chronicle and graph views, and the disaster triggers.

## 11. Validation

Automated behaviour tests of the promises: a faction on barren land dies out;
population follows food capacity; scarcity raises the chance of war; the same
seed and log give the same history.

Plus one look at Earth: the food map computed from Earth's climate should put
people where early civilisations were (the Nile, the Ganges, the Yellow River).
That is for tuning by eye, never a pass/fail test, and there is no comparison
with real historical statistics.

## 12. Stages

Building starts no earlier than `fms-editor`'s first milestone (the shared
render and editing code) and `fms-climate` M3 step 0 (the `.fmsmap` crate).
Headless, with scripted setups and PNG dumps, until the GUI stage.

| Stage | Contents |
|---|---|
| Provinces and food | Province generation, food capacity with its fallbacks |
| Population | Growth, famine, settling, the sea |
| War | Wars, capitals, conquest |
| Culture and religion | Shares, drift, breakaways |
| Disasters | Natural rates from geo and climate layers, effects |
| Chronicle and names | Events, name generator, graphs data |
| GUI | Drawing, trait panels, timeline, disaster triggers, mid-run edits |
| Export | `.fmsmap` layers and name tables, project file |

## 13. Reserved for later

- **Technology and ages.**
- **Trade.**
- **Diplomacy and alliances.**
- **Rulers and dynasties.**
- **Armies and navies.**
- **An economy beyond food**: materials such as metal or salt would need
  `fms-geo`'s mineral layer.
- **Plague.**
- **Land use and its effect on climate.**
- **Branching timelines.**
- **A world pre-filled with stateless peoples.**
- **New cultures and religions arising by themselves.**
