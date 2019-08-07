------------------------------------------------------------
-- Laden von benötigten Skripten
------------------------------------------------------------
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

------------------------------------------------------------
-- Parsen der Parameter
------------------------------------------------------------

-- Verwendetes Modell
geometrieName	= util.GetParam("-grid", "Modell1", "Dateiname der Geometrie")

-- Anzahl der Verfeinerungen
verfeinerungen		= util.GetParamNumber("-numRefs", 1, "Anzahl der Verfeinerungen")

-- Diffusionstensor, nur erforderlich für Modell2
diffusion		= util.GetParamNumber("-diff", 0.1, "Werte auf der Hauptdiagonale")

------------------------------------------------------------
-- Laden des Gebiets und dessen Verfeinerungen
------------------------------------------------------------

function ErstelleGebiet() 
	InitUG(3, AlgebraType("CPU", 1));

	-- Definition der benötigten Subsets
	if  geometrieName == "Modell1" then
		benoetigteSubsets = {"innen", "oben", "unten"}
	elseif  geometrieName == "Modell2" then
		benoetigteSubsets = {"innen", "membran","oben", "unten"}
	else
		print("Es gibt leider keine passende Geometrie mit diesem Namen. Bitte geben Sie entweder Modell1, Modell2 oder nichts als Parameter für -grid ein.")
	end

	-- Laden des Gebiets
	gebiet = util.CreateDomain("grids/"..geometrieName..".ugx", 0, benoetigteSubsets)

	-- Verfeinerungen des Gebiets
	util.refinement.CreateRegularHierarchy(gebiet, verfeinerungen, true)

	return gebiet

end

------------------------------------------------------------
-- Erstellung des Approximationsraums
------------------------------------------------------------

function ErstelleApproximationsraum(gebiet)

	approxRaum = ApproximationSpace(gebiet)
	approxRaum:add_fct("c", "Lagrange", 1)

	approxRaum:init_levels()
	approxRaum:init_top_surface()

	return approxRaum
end

------------------------------------------------------------
-- Diskretisierung
------------------------------------------------------------
function DiffTensorMembran()
	return	diffusion, 0, 0,
			0, diffusion, 0,
			0, 0, diffusion
end

function ErstelleDiskretisierung(approxRaum)
	elemDisk = {}
	elemDisk["Innen"] = ConvectionDiffusion("c", "innen", "fv1")
	elemDisk["Innen"]:set_diffusion(1.0)
	elemDisk["Innen"]:set_source(0)

	if  geometrieName == "Modell2" then
		elemDisk["Membran"] = ConvectionDiffusion("c", "membran", "fv1")
		elemDisk["Membran"]:set_diffusion("DiffTensorMembran")
		elemDisk["Membran"]:set_source(0)
	end

	-- Dirichlet-Randbedingung
	dirichletRB = DirichletBoundary()
	dirichletRB:add(1.0, "c", "oben")

	-- Neumann-Randbedingung
	neumannRB = NeumannBoundary("c")
	neumannRB:add(1.0, "unten", "innen")


	disk = DomainDiscretization(approxRaum)
	disk:add(elemDisk["Innen"])
	if  geometrieName == "Modell2" then
		disk:add(elemDisk["Membran"])
	end
	disk:add(dirichletRB)
	disk:add(neumannRB)

	return disk
end

------------------------------------------------------------
-- Erstellung des Lösers
------------------------------------------------------------
function ErstelleLoeser(approxRaum)

	loeser = util.solver.CreateSolver({
		type = "bicgstab",
		precond = {
			type		= "gmg",
			approxSpace	= approxRaum,
			smoother	= "jac",
			baseSolver	= "lu"
		}})

	return loeser
end 

------------------------------------------------------------
-- Hauptprogramm
------------------------------------------------------------

print("Skript simulationDiffusion wird gestartet")


local gebiet = ErstelleGebiet()
local approxRaum = ErstelleApproximationsraum(gebiet)
local disk = ErstelleDiskretisierung(approxRaum)
local loeser = ErstelleLoeser(approxRaum)

A = AssembledLinearOperator(disk)
u = GridFunction(approxRaum)
b = GridFunction(approxRaum)
u:set(0.0)

disk:adjust_solution(u)
disk:assemble_linear(A, b)

loeser:init(A, u)
loeser:apply(u, b)


dateiName = "Verfeinerung" .. verfeinerungen .. geometrieName

if geometrieName == "Modell2" then
	dateiName = dateiName .. "Diffusion" .. string.sub(diffusion, 3)
end
print("schreibe Lösung nach " .. dateiName)
WriteGridFunctionToVTK(u, dateiName)


print("Skript simulationDiffusion ist beendet")
