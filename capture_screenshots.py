import sys
import os
import time
from pathlib import Path
import pandas as pd
from PySide6.QtWidgets import QApplication, QHeaderView
from PySide6.QtCore import Qt, QTimer

# Add project root to sys.path
root_dir = Path(r'E:\Github\SAM_Header')
sys.path.insert(0, str(root_dir))

from GUI.app import MainWindow

# Set environment variables for offscreen rendering
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
os.environ['QT_QPA_FONTDIR'] = r'C:\Windows\Fonts'

def create_demo_csv(path, is_block=False):
    data = [
        ['site_info', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['site_id', 'Test Site', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['site_name', 'Demo Site', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['site_lat', '45.0', '(ºN)', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['site_long', '-120.0', '(ºE)', '', '', '', '', 'Sun Compass Information', '', '', '', '', '', '', '', '', 'Calculated Fields', '', '', '', '', ''],
        ['site_elevation', '100', '(meters)', 'only used if no sun data', '', '', '', '[default is yes] (yes or no)', 'all sun compass info is optional', '', '', '', '', '', '', '', '', 'Optional Field', 'default core strike', '', '', ''],
        ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['sample_name', 'comment', 'strat_level', 'magnetic_core_strike', 'core_dip', 'bedding_strike', 'bedding_dip', 'correct_bedding_using_local_dec', 'shadow_angle', 'GMT_offset', 'year', 'month', 'days', 'hours', 'minutes', 'mass', 'runs', 'sun_core_strike', 'calculated_IGRF', 'IGRF_local_dec', 'calculated_mag_dec', 'core_strike', 'corrected_bedding_strike']
    ]
    
    if is_block:
        # For block samples, use negative shadow angle as per instructions
        data.append(['B1', 'test block', '10', '120', '0', '30', '45', 'yes', '-150', '-8', '2023', '05', '20', '10', '30', '1', '1', '', '', '', '', '', ''])
    else:
        data.append(['C1', 'test core', '5', '45', '10', '30', '45', 'yes', '150', '-8', '2023', '05', '20', '10', '30', '1', '1', '', '', '', '', '', ''])
        data.append(['C2', 'test core 2', '6', '180', '15', '30', '45', 'yes', '160', '-8', '2023', '05', '20', '11', '00', '1', '1', '', '', '', '', '', ''])
    
    df = pd.DataFrame(data)
    df.to_csv(path, index=False, header=False)
    return path

def capture_startup(app, screenshot_dir):
    window = MainWindow()
    window.show()
    app.processEvents()
    
    # Startup capture
    path = screenshot_dir / 'startup.png'
    window.grab().save(str(path))
    size = window.size()
    window.close()
    return path, size

def capture_loaded(app, screenshot_dir, csv_path, filename, is_block=False):
    window = MainWindow()
    window.show()
    app.processEvents()
    
    # Simulate adding files
    window._add_files([str(csv_path)])
    
    if is_block:
        # Find the orientation mode combo box and set it to block
        # Based on GUI code exploration, we might need to find the widget
        for combo in window.findChildren(objectName='orientation_mode_combo'): # Check if it has a name
             combo.setCurrentText('block')
        # Or wait, let's look for how it's handled in the UI
        # By default it's likely pomeroy/core
        pass

    app.processEvents()
    
    # Resize window to fit table
    # Find the table widget in the preview panel
    table = None
    for widget in window.findChildren(objectName='orientation_table'):
        table = widget
        break
    
    if table:
        # Calculate width
        width = table.verticalHeader().width() + 40 # margin
        for i in range(table.columnCount()):
            width += table.columnWidth(i)
        
        # Add width for other panels
        # The preview panel is 5/11 of the total width roughly
        current_size = window.size()
        new_width = int(width * (11/5) + 100) # conservative estimate
        window.resize(new_width, current_size.height())
        app.processEvents()

    path = screenshot_dir / filename
    window.grab().save(str(path))
    size = window.size()
    window.close()
    return path, size

if __name__ == '__main__':
    app = QApplication(sys.argv)
    screenshot_dir = root_dir / 'docs' / 'screenshots'
    screenshot_dir.mkdir(parents=True, exist_ok=True)
    
    core_csv = create_demo_csv(root_dir / 'demo_core.csv', is_block=False)
    block_csv = create_demo_csv(root_dir / 'demo_block.csv', is_block=True)
    
    results = []
    
    # Startup
    p1, s1 = capture_startup(app, screenshot_dir)
    results.append(f'startup.png: {s1.width()}x{s1.height()}')
    
    # Loaded Core
    p2, s2 = capture_loaded(app, screenshot_dir, core_csv, 'loaded-core.png', is_block=False)
    results.append(f'loaded-core.png: {s2.width()}x{s2.height()}')
    
    # Loaded Block
    # For block, we need to ensure the mode is switched.
    # Looking at app.py _build_options_panel, there is a combo for mode.
    p3, s3 = capture_loaded(app, screenshot_dir, block_csv, 'loaded-block.png', is_block=True)
    results.append(f'loaded-block.png: {s3.width()}x{s3.height()}')
    
    for r in results:
        print(r)
    
    # Cleanup
    os.remove(core_csv)
    os.remove(block_csv)

