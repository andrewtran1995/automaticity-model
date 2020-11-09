function saveOutput(config,neurons)
    save(sprintf("output\\%s-%s", class(config), datestr(now, "yyyy-mm-dd-HH-MM")), "config", "neurons");
end

