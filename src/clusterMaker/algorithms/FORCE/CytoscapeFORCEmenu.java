package clusterMaker.algorithms.FORCE;

import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import cytoscape.Cytoscape;
import cytoscape.view.CytoscapeDesktop;
import cytoscape.view.cytopanels.CytoPanel;


public class CytoscapeFORCEmenu extends JMenu implements ActionListener {
	
	private CytoscapeFORCEgui gui = null;
	private JPanel guiDummyPanel = null;
	
	private JMenuItem startGuiButton;
	private JMenuItem startGuiPanelButton;
	private JMenuItem stopGuiButton;
	
	private JScrollPane guiDummyPanelScrollPane;
	
	public CytoscapeFORCEmenu() {
		super("FORCE clustering");
		
		addMenuItems();
	}
	
	private void addMenuItems() {
		
		startGuiButton = new JMenuItem("Start FORCE clustering GUI");
		startGuiButton.setActionCommand("embedd");
		startGuiButton.addActionListener(this);
		
		startGuiPanelButton = new JMenuItem("Start FORCE GUI in an own frame");
		startGuiPanelButton.setActionCommand("frame");
		startGuiPanelButton.addActionListener(this);
		
		stopGuiButton = new JMenuItem("Remove/stop FORCE clustering GUI");
		stopGuiButton.setActionCommand("remove");
		stopGuiButton.addActionListener(this);
		stopGuiButton.setEnabled(false);
		
		JMenuItem infoButton = new JMenuItem("About/help");
		infoButton.setActionCommand("info");
		infoButton.addActionListener(this);
		
		this.add(startGuiButton);
		this.add(startGuiPanelButton);
		this.add(stopGuiButton);
		this.add(infoButton);
		
	}
	
	public void actionPerformed(ActionEvent e) {
		
		String c = e.getActionCommand();
		
		if (c.equalsIgnoreCase("embedd")) {
			
			embeddGui();
			
		} else if (c.equalsIgnoreCase("frame")) {
			
			startGuiFrame();
			
		} else if (c.equalsIgnoreCase("remove")) {
			
			removeGui();
			
		} else if (c.equalsIgnoreCase("info")) {
			
			new CytoscapeFORCEinfoFrame();
			
		}
		
	}
	
	private void startGuiFrame() {
		
		gui = new CytoscapeFORCEgui();
		
		JFrame f = new JFrame("FORCE");
		f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		f.getContentPane().add(new JScrollPane(gui));
		
		f.setVisible(true);
		f.pack();
		
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		int x = (screenSize.width-f.getSize().width)/2;
		int y = (screenSize.height-f.getSize().height)/2;
		f.setLocation(x, y);
	}
	
	private void embeddGui() {
		
		if (this.gui == null) {
			CytoscapeDesktop desktop = Cytoscape.getDesktop();
			CytoPanel cytoPanel = desktop.getCytoPanel (SwingConstants.WEST);
			gui = new CytoscapeFORCEgui();
			
			guiDummyPanel = new JPanel();
			guiDummyPanel.add(gui);
			guiDummyPanelScrollPane = new JScrollPane(guiDummyPanel);
			
			cytoPanel.add("FORCE Clustering", guiDummyPanelScrollPane);
			
			startGuiButton.setEnabled(false);
			stopGuiButton.setEnabled(true);
			
		}
		
	}
	
	private void removeGui() {
		
		if (this.gui != null) {
			CytoscapeDesktop desktop = Cytoscape.getDesktop();
			CytoPanel cytoPanel = desktop.getCytoPanel (SwingConstants.WEST);
			
			cytoPanel.remove(guiDummyPanelScrollPane);
			gui = null;
			
			startGuiButton.setEnabled(true);
			stopGuiButton.setEnabled(false);
			
		}
		
	}
	
	
}