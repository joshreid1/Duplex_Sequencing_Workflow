function split(str, sep) -- splits string by delimiter
	local sep, fields = sep or ":", {}
	local pattern = string.format("([^%s]+)", sep)
	str:gsub(pattern, function(c) fields[#fields+1] = c end)
	return fields
end

function indexOf(table, value) -- returns the first index of the table that matched the value
	for i, v in ipairs(table) do
		if v == value then
			return i
		end
	end
	return nil
end

function spliceai(entry) -- processes SpliceAI scores to check for maximum transcript score
	local t = type(entry)
	if t == "string" then
		return entry -- returns original value if single entry
	elseif t == "table" then
		local maximums = {}
		for i=1,#entry do -- calculate the maximum SpliceAI score of AG, AL, DG, DL for each entry
			maximums[i] = math.max(tonumber(split(entry[i], "|")[3]), tonumber(split(entry[i], "|")[4]), tonumber(split(entry[i], "|")[5]), tonumber(split(entry[i], "|")[6]))
		end
		return entry[indexOf(maximums, math.max(unpack(maximums)))] -- returns the full record which contains the maximum SpliceAI score
	end
end
