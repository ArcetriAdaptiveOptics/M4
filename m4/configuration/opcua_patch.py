import re
from enum import IntEnum
from datetime import datetime
import uuid
from opcua import ua

try:
    from opcua.common import structures


    def patched_generate_python_class(model, env=None):
        if env is None: env = {}
        if "ua" not in env: env['ua'] = ua
        if "datetime" not in env: env['datetime'] = datetime
        if "uuid" not in env: env['uuid'] = uuid
        if "IntEnum" not in env: env['IntEnum'] = IntEnum

        for element in model:
            code = element.get_code()
            code = re.sub(r"(\s*)(\w+[/\-]\w+)(\s*=\s*\d+)",
                          lambda m: m.group(1) + m.group(2).replace('/', '_').replace('-', '_') + m.group(3), code)

            # The 'exec' is kept, but error printing is removed.
            # The 'raise' will still stop the program on a critical error.
            exec(code, env)

        return env


    structures._generate_python_class = patched_generate_python_class

except ImportError:
    # Errors are no longer printed to the console.
    pass
