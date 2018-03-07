from chemProfileLib.cli import t
from chemProfileLib.tools import tool_list
from chemProfileLib.tools.exceptions import ToolNotFoundException

def check():
    """
    Checks if all tools are installed and if needed indecies (genome) are available.
    :return:
    """
    print("Checking if all tools are installed.")
    for tool in tool_list:
        try:
            tool_list[tool].check()
            state = True
        except ToolNotFoundException:
            state = False

        print(" - {:<20} {:<15}".format(tool, t.green("Installed") if state is True else t.red("Not found")))